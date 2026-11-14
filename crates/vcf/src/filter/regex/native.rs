//! C-locale REG_NOSUB semantics for patterns containing back references.

use std::ffi::CString;
use std::fmt;
use std::mem::MaybeUninit;
use std::sync::{Arc, Mutex};

pub(super) const LIMIT: usize = 1_000_000;

struct Locale(libc::locale_t);

impl Locale {
    fn new() -> Option<Self> {
        // SAFETY: constant mask/name and a null base request a new locale.
        // The returned non-null handle is owned until this object's Drop.
        let raw =
            unsafe { libc::newlocale(libc::LC_ALL_MASK, c"C".as_ptr(), std::ptr::null_mut()) };
        if raw.is_null() { None } else { Some(Self(raw)) }
    }
}

impl Drop for Locale {
    fn drop(&mut self) {
        // SAFETY: this handle came from successful newlocale, is not the
        // global-locale sentinel, and has no outstanding LocaleScope borrow.
        unsafe { libc::freelocale(self.0) }
    }
}

struct LocaleScope<'a> {
    previous: libc::locale_t,
    _locale: &'a Locale,
}

impl<'a> LocaleScope<'a> {
    fn enter(locale: &'a Locale) -> Option<Self> {
        // SAFETY: locale is a live owned handle. uselocale changes only this
        // thread; the borrow prevents freeing it before scope restoration.
        let previous = unsafe { libc::uselocale(locale.0) };
        if previous.is_null() {
            None
        } else {
            Some(Self {
                previous,
                _locale: locale,
            })
        }
    }
}

impl Drop for LocaleScope<'_> {
    fn drop(&mut self) {
        // SAFETY: previous is the caller's still-live locale or the valid
        // LC_GLOBAL_LOCALE sentinel. It is restored, never freed here.
        // The raw-pointer fields make this private guard non-Send.
        unsafe {
            libc::uselocale(self.previous);
        }
    }
}

struct OwnedRegex {
    raw: Box<libc::regex_t>,
    locale: Locale,
}

// SAFETY: POSIX compiled regex/locale handles have no thread affinity. The
// object stays at its stable Box address. Every execution and final cleanup
// is exclusive under its owning Mutex; no raw handle escapes. LocaleScope
// stays on the calling thread and restores its previous locale before unlock.
// No unsafe Sync is needed: Mutex supplies shared access from this Send bound.
unsafe impl Send for OwnedRegex {}

impl Drop for OwnedRegex {
    fn drop(&mut self) {
        let _scope = LocaleScope::enter(&self.locale);
        // SAFETY: regcomp succeeded exactly once at this stable address;
        // final ownership precludes concurrent use and regfree runs once.
        // regfree only releases the compiled allocation; if locale activation
        // fails, cleanup itself has no locale-dependent matching operation.
        unsafe { libc::regfree(&mut *self.raw) }
        // Scope restores the thread locale before Locale's field destructor.
    }
}

#[derive(Clone)]
pub(super) struct NativeRegex {
    owned: Arc<Mutex<OwnedRegex>>,
    program_size: usize,
}

impl fmt::Debug for NativeRegex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("NativeRegex")
            .field("program_size", &self.program_size)
            .finish_non_exhaustive()
    }
}

fn resource(what: &'static str) {
    super::super::undefined(what);
}

impl NativeRegex {
    pub(super) fn new(
        pattern: &[u8],
        icase: bool,
        program_size: usize,
    ) -> Result<Self, &'static str> {
        if program_size == 0 || program_size > LIMIT || pattern.len() > LIMIT {
            resource("native regex expanded-program size limit exceeded");
            return Err("REG_ESPACE");
        }
        let pattern = CString::new(pattern).map_err(|_| "REG_BADPAT")?;
        let Some(locale) = Locale::new() else {
            resource("native regex C locale allocation failed");
            return Err("REG_ESPACE");
        };
        // Keep the C object at one heap address from compilation to regfree.
        let mut storage = Box::new(MaybeUninit::<libc::regex_t>::zeroed());
        let raw = {
            let Some(_scope) = LocaleScope::enter(&locale) else {
                resource("native regex C locale activation failed");
                return Err("REG_ESPACE");
            };
            let flags = libc::REG_NOSUB | if icase { libc::REG_ICASE } else { 0 };
            // SAFETY: storage is aligned, sized, stable writable regex_t
            // storage; pattern is a live terminated C string. No other
            // operation can observe this unpublished object.
            let result = unsafe { libc::regcomp(storage.as_mut_ptr(), pattern.as_ptr(), flags) };
            if result != 0 {
                if result == libc::REG_ESPACE {
                    resource("native regex compilation ran out of space");
                }
                // A failed regcomp object is not passed to regfree.
                return Err("REG_BADPAT");
            }
            // SAFETY: regcomp returned success and initialized regex_t;
            // assume_init moves only the Box handle, not its allocation.
            unsafe { storage.assume_init() }
        };
        Ok(Self {
            owned: Arc::new(Mutex::new(OwnedRegex { raw, locale })),
            program_size,
        })
    }

    pub(super) fn is_match(&self, subject: &[u8]) -> bool {
        let subject = &subject[..subject
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(subject.len())];
        if subject
            .len()
            .checked_add(1)
            .and_then(|n| n.checked_mul(self.program_size))
            .filter(|&n| n <= LIMIT)
            .is_none()
        {
            resource("native regex subject/program size limit exceeded");
            return false;
        }
        let Ok(subject) = CString::new(subject) else {
            resource("native regex subject is not a C string");
            return false;
        };
        let Ok(owned) = self.owned.lock() else {
            resource("native regex mutex poisoned");
            return false;
        };
        let Some(_scope) = LocaleScope::enter(&owned.locale) else {
            resource("native regex C locale activation failed");
            return false;
        };
        // SAFETY: the compiled object is exclusively locked and live; both
        // compilation and matching use its owned C locale. REG_NOSUB with
        // nmatch=0 requires no output array. The input is NUL terminated.
        match unsafe { libc::regexec(&*owned.raw, subject.as_ptr(), 0, std::ptr::null_mut(), 0) } {
            0 => true,
            libc::REG_NOMATCH => false,
            _ => {
                resource("native regex execution failed");
                false
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::NativeRegex;

    fn current_locale() -> libc::locale_t {
        // SAFETY: a null argument queries this thread's locale without
        // modifying it. The returned caller-owned handle is never freed.
        unsafe { libc::uselocale(std::ptr::null_mut()) }
    }

    #[test]
    fn failed_compilation_size_rejection_and_shared_drop_restore_state() {
        super::super::super::take_undefined();
        let before = current_locale();
        assert!(NativeRegex::new(b"[", false, 1).is_err());
        assert_eq!(current_locale(), before);
        let native = NativeRegex::new(b"^\\(a\\)\\1$", false, 16).unwrap();
        let oversized = b"a".repeat(65_537);
        assert!(!native.is_match(&oversized));
        assert!(super::super::super::take_undefined().is_some());
        assert_eq!(current_locale(), before);
        assert!(native.is_match(b"aa\0ignored"));
        assert!(!native.is_match(b"ab\0aa"));
        assert!(super::super::super::take_undefined().is_none());
        let shared = native.clone();
        drop(native);
        assert!(shared.is_match(b"aa"));
        drop(shared);
        assert_eq!(current_locale(), before);
    }

    #[test]
    fn poisoned_native_handle_reports_unsupported_instead_of_a_match() {
        super::super::super::take_undefined();
        let before = current_locale();
        let native = NativeRegex::new(b"^\\(a\\)\\1$", false, 16).unwrap();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _lock = native.owned.lock().unwrap();
            panic!("poison the native regex lock");
        }));
        assert!(result.is_err());
        assert!(!native.is_match(b"aa"));
        assert!(super::super::super::take_undefined().is_some());
        drop(native);
        assert_eq!(current_locale(), before);
    }
}
