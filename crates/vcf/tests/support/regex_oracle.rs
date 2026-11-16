use std::ffi::CString;
use std::mem::MaybeUninit;

pub fn matches(pattern: &str, subject: &str, icase: bool) -> Result<bool, ()> {
    let pattern = CString::new(pattern).unwrap();
    let subject = CString::new(subject).unwrap();
    let mut regex = MaybeUninit::<libc::regex_t>::uninit();
    let flags = libc::REG_NOSUB | if icase { libc::REG_ICASE } else { 0 };
    // SAFETY: regcomp initializes regex_t on success, and only that branch
    // passes it to regexec/regfree. Both inputs are live terminated C strings;
    // REG_NOSUB and zero nmatch require no output match array. Nothing borrows
    // the C allocation while regfree releases it.
    unsafe {
        if libc::regcomp(regex.as_mut_ptr(), pattern.as_ptr(), flags) != 0 {
            return Err(());
        }
        let mut regex = regex.assume_init();
        let status = libc::regexec(&regex, subject.as_ptr(), 0, std::ptr::null_mut(), 0);
        libc::regfree(&mut regex);
        match status {
            0 => Ok(true),
            libc::REG_NOMATCH => Ok(false),
            _ => panic!("regexec failed: {status}"),
        }
    }
}
