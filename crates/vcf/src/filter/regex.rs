//! POSIX basic regular expressions as glibc's `regcomp(REG_NOSUB [|
//! REG_ICASE])` and `regexec` read and match them in the C locale: the
//! syntax of `RE_SYNTAX_POSIX_BASIC` with the GNU operators (`\+ \? \|
//! \w \W \s \S \b \B \< \> \` \'`), back-references, bracket expressions
//! with classes, collating symbols and equivalence classes, and the same
//! compile errors. Only whether a match exists is reported.
//!
//! Under `REG_ICASE` glibc upper-cases the pattern and the subject and
//! compares them byte for byte; characters taken after a backslash and
//! class names are read from the pattern as written.
//!
//! Short, branch-free patterns use the Rust matcher even with back references.
//! Other back-reference patterns use system REG_NOSUB in a private C locale.

use std::cell::Cell;

mod native;
use native::NativeRegex;

const DUP_MAX: i64 = 0x7fff;
const TREE_DEPTH: usize = 64;
const MATCH_DEPTH: usize = 128;
const MATCH_WORK: usize = 1_000_000;

#[derive(Clone, Copy)]
struct Budget {
    work: usize,
    depth: usize,
    exhausted: bool,
}

thread_local! {
    static BUDGET: Cell<Budget> = const { Cell::new(Budget { work: 0, depth: 0, exhausted: false }) };
}

fn charge(work: usize) -> bool {
    BUDGET.with(|state| {
        let mut b = state.get();
        if b.exhausted || work > b.work {
            b.exhausted = true;
            state.set(b);
            super::undefined("regular expression work budget exceeded");
            return false;
        }
        b.work -= work;
        state.set(b);
        true
    })
}

fn exhausted() -> bool {
    BUDGET.with(|state| state.get().exhausted)
}

struct Frame;
impl Frame {
    fn enter() -> Option<Self> {
        if !charge(1) {
            return None;
        }
        BUDGET.with(|state| {
            let mut b = state.get();
            if b.depth == MATCH_DEPTH {
                b.exhausted = true;
                state.set(b);
                super::undefined("regular expression recursion budget exceeded");
                return None;
            }
            b.depth += 1;
            state.set(b);
            Some(Frame)
        })
    }
}
impl Drop for Frame {
    fn drop(&mut self) {
        BUDGET.with(|state| {
            let mut b = state.get();
            b.depth -= 1;
            state.set(b);
        });
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Anchor {
    LineFirst,
    LineLast,
    BufFirst,
    BufLast,
    WordFirst,
    WordLast,
    WordDelim,
    NotWordDelim,
}

#[derive(Debug, Clone)]
enum Node {
    Empty,
    Char(u8),
    Any,
    Set(Box<[bool; 256]>),
    Assert(Anchor),
    Group(Box<Node>, usize),
    BackRef(usize),
    Concat(Vec<Node>),
    Alt(Vec<Node>),
    Repeat(Box<Node>, u32, Option<u32>),
}

impl Node {
    fn depth(&self) -> usize {
        match self {
            Node::Group(n, _) | Node::Repeat(n, _, _) => n.depth() + 1,
            Node::Concat(v) | Node::Alt(v) => v.iter().map(Node::depth).max().unwrap_or(0) + 1,
            _ => 1,
        }
    }

    /// Children are checked before attaching their parent, so inspecting and
    /// dropping this node visits at most TREE_DEPTH + 1 recursive levels.
    fn checked(self) -> Result<Self, &'static str> {
        if self.depth() > TREE_DEPTH {
            super::undefined("regular expression syntax tree depth budget exceeded");
            Err("REG_ESPACE")
        } else {
            Ok(self)
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Tok {
    Char(u8),
    BackSlash,
    Alt,
    BackRef(usize),
    Anchor(Anchor),
    Word(bool),
    Space(bool),
    OpenSub,
    CloseSub,
    DupAsterisk,
    DupPlus,
    DupQuestion,
    OpenDup,
    CloseDup,
    OpenBracket,
    Period,
    End,
}

/// A compiled pattern.
#[derive(Debug, Clone)]
pub(crate) struct Regex {
    root: Node,
    icase: bool,
    groups: usize,
    literal: Option<Literal>,
    native: Option<NativeRegex>,
}

#[derive(Debug, Clone)]
struct Literal {
    bytes: Vec<u8>,
    failure: Vec<usize>,
    first: bool,
    last: bool,
}

impl Literal {
    fn new(root: &Node) -> Option<Self> {
        let nodes = match root {
            Node::Concat(v) => &v[..],
            _ => std::slice::from_ref(root),
        };
        let mut bytes = Vec::new();
        let (mut first, mut last) = (false, false);
        for (i, node) in nodes.iter().enumerate() {
            match node {
                Node::Char(c) => bytes.push(*c),
                Node::Assert(Anchor::LineFirst | Anchor::BufFirst) if i == 0 => first = true,
                Node::Assert(Anchor::LineLast | Anchor::BufLast) if i + 1 == nodes.len() => {
                    last = true
                }
                _ => return None,
            }
        }
        let mut failure = vec![0; bytes.len()];
        let mut j = 0;
        for i in 1..bytes.len() {
            while j > 0 && bytes[i] != bytes[j] {
                j = failure[j - 1];
            }
            if bytes[i] == bytes[j] {
                j += 1;
            }
            failure[i] = j;
        }
        Some(Self {
            bytes,
            failure,
            first,
            last,
        })
    }

    fn is_match(&self, re: &Regex, s: &[u8]) -> bool {
        let n = self.bytes.len();
        let same = |slice: &[u8]| slice.iter().zip(&self.bytes).all(|(&c, &p)| re.ch(c) == p);
        if self.first && self.last {
            return s.len() == n && same(s);
        }
        if self.first {
            return s.len() >= n && same(&s[..n]);
        }
        if self.last {
            return s.len() >= n && same(&s[s.len() - n..]);
        }
        if n == 0 {
            return true;
        }
        let mut j = 0;
        for &c in s {
            let c = re.ch(c);
            while j > 0 && c != self.bytes[j] {
                j = self.failure[j - 1];
            }
            if c == self.bytes[j] {
                j += 1;
            }
            if j == n {
                return true;
            }
        }
        false
    }
}

struct Parser<'a> {
    /// The pattern as glibc's parser sees it (upper-cased under ICASE).
    buf: Vec<u8>,
    /// The pattern as written.
    raw: &'a [u8],
    at: usize,
    tok: Tok,
    tok_len: usize,
    groups: usize,
    completed: u32,
    icase: bool,
}

fn is_word(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_'
}

impl Parser<'_> {
    fn byte(&self, i: usize) -> Option<u8> {
        self.buf.get(i).copied()
    }

    /// `peek_token` at the current position.
    fn peek(&self, caret_here: bool) -> (Tok, usize) {
        self.peek_at(self.at, caret_here)
    }

    fn peek_at(&self, at: usize, caret_here: bool) -> (Tok, usize) {
        let Some(c) = self.byte(at) else {
            return (Tok::End, 0);
        };
        if c == b'\\' {
            let Some(&c2) = self.raw.get(at + 1) else {
                return (Tok::BackSlash, 1);
            };
            let t = match c2 {
                b'|' => Tok::Alt,
                b'1'..=b'9' => Tok::BackRef((c2 - b'1') as usize),
                b'<' => Tok::Anchor(Anchor::WordFirst),
                b'>' => Tok::Anchor(Anchor::WordLast),
                b'b' => Tok::Anchor(Anchor::WordDelim),
                b'B' => Tok::Anchor(Anchor::NotWordDelim),
                b'w' => Tok::Word(false),
                b'W' => Tok::Word(true),
                b's' => Tok::Space(false),
                b'S' => Tok::Space(true),
                b'`' => Tok::Anchor(Anchor::BufFirst),
                b'\'' => Tok::Anchor(Anchor::BufLast),
                b'(' => Tok::OpenSub,
                b')' => Tok::CloseSub,
                b'+' => Tok::DupPlus,
                b'?' => Tok::DupQuestion,
                b'{' => Tok::OpenDup,
                b'}' => Tok::CloseDup,
                _ => Tok::Char(c2),
            };
            return (t, 2);
        }
        let t = match c {
            b'*' => Tok::DupAsterisk,
            b'[' => Tok::OpenBracket,
            b'.' => Tok::Period,
            b'^' if at == 0 || caret_here => Tok::Anchor(Anchor::LineFirst),
            b'$' => {
                if at + 1 == self.buf.len() {
                    Tok::Anchor(Anchor::LineLast)
                } else {
                    // A following '$' cannot be an alternation/group-close;
                    // avoid recursive lookahead over a run of dollar signs.
                    let next = if self.byte(at + 1) == Some(b'$') {
                        Tok::Char(b'$')
                    } else {
                        self.peek_at(at + 1, false).0
                    };
                    if matches!(next, Tok::Alt | Tok::CloseSub) {
                        Tok::Anchor(Anchor::LineLast)
                    } else {
                        Tok::Char(c)
                    }
                }
            }
            _ => Tok::Char(c),
        };
        (t, 1)
    }

    /// `fetch_token`: consumes the current token and reads the next.
    fn fetch(&mut self, caret_here: bool) {
        self.at += self.tok_len;
        let (t, n) = self.peek(caret_here);
        self.tok = t;
        self.tok_len = n;
    }

    fn reg_exp(&mut self, nest: usize) -> Result<Node, &'static str> {
        let initial = self.completed;
        let mut alts = vec![self.branch(nest)?];
        while self.tok == Tok::Alt {
            self.fetch(true);
            if self.tok != Tok::Alt
                && self.tok != Tok::End
                && (nest == 0 || self.tok != Tok::CloseSub)
            {
                let accumulated = self.completed;
                self.completed = initial;
                alts.push(self.branch(nest)?);
                self.completed |= accumulated;
            } else {
                alts.push(Node::Empty);
            }
        }
        (if alts.len() == 1 {
            alts.pop().unwrap()
        } else {
            Node::Alt(alts)
        })
        .checked()
    }

    fn branch(&mut self, nest: usize) -> Result<Node, &'static str> {
        let mut seq = vec![self.expression(nest)?];
        while self.tok != Tok::Alt
            && self.tok != Tok::End
            && (nest == 0 || self.tok != Tok::CloseSub)
        {
            seq.push(self.expression(nest)?);
        }
        (if seq.len() == 1 {
            seq.pop().unwrap()
        } else {
            Node::Concat(seq)
        })
        .checked()
    }

    fn expression(&mut self, nest: usize) -> Result<Node, &'static str> {
        let mut tree = match self.tok {
            Tok::Char(c) => Node::Char(c),
            Tok::OpenSub => self.sub_exp(nest + 1)?,
            Tok::OpenBracket => self.bracket()?,
            Tok::BackRef(i) => {
                if self.completed & (1 << i) == 0 {
                    return Err("REG_ESUBREG");
                }
                Node::BackRef(i)
            }
            Tok::OpenDup => return Err("REG_BADRPT"),
            // At the start of an expression these are ordinary.
            Tok::DupAsterisk => Node::Char(b'*'),
            Tok::DupPlus => Node::Char(b'+'),
            Tok::DupQuestion => Node::Char(b'?'),
            Tok::CloseSub => return Err("REG_ERPAREN"),
            Tok::CloseDup => Node::Char(b'}'),
            Tok::Anchor(a) => {
                // Anchors take no repetition: what follows starts afresh.
                self.fetch(false);
                return Ok(Node::Assert(a));
            }
            Tok::Period => Node::Any,
            Tok::Word(neg) => Node::Set(Box::new(class_set(is_word, neg))),
            Tok::Space(neg) => Node::Set(Box::new(class_set(crate::cnum::is_space, neg))),
            Tok::Alt | Tok::End => return Ok(Node::Empty),
            Tok::BackSlash => return Err("REG_EESCAPE"),
        };
        self.fetch(false);
        while matches!(
            self.tok,
            Tok::DupAsterisk | Tok::DupPlus | Tok::DupQuestion | Tok::OpenDup
        ) {
            tree = self.dup(tree)?;
            tree = tree.checked()?;
            // In a BRE no `*` or `\{` may follow a repetition.
            if matches!(self.tok, Tok::DupAsterisk | Tok::OpenDup) {
                return Err("REG_BADRPT");
            }
        }
        Ok(tree)
    }

    fn sub_exp(&mut self, nest: usize) -> Result<Node, &'static str> {
        if nest > TREE_DEPTH {
            super::undefined("regular expression parser depth budget exceeded");
            return Err("REG_ESPACE");
        }
        let index = self.groups;
        self.groups += 1;
        self.fetch(true);
        let inner = if self.tok == Tok::CloseSub {
            Node::Empty
        } else {
            let inner = self.reg_exp(nest)?;
            if self.tok != Tok::CloseSub {
                return Err("REG_EPAREN");
            }
            inner
        };
        if index <= 8 {
            self.completed |= 1 << index;
        }
        Node::Group(Box::new(inner), index).checked()
    }

    /// `fetch_number`: -1 for none, -2 for invalid or the end.
    fn number(&mut self) -> i64 {
        let mut num = -1i64;
        loop {
            self.fetch(false);
            let c = match self.tok {
                Tok::Char(c) => Some(c),
                _ => None,
            };
            if self.tok == Tok::End {
                return -2;
            }
            if self.tok == Tok::CloseDup || self.tok_char() == Some(b',') {
                break;
            }
            num = match c {
                Some(d @ b'0'..=b'9') if num != -2 => {
                    if num == -1 {
                        i64::from(d - b'0')
                    } else {
                        (num * 10 + i64::from(d - b'0')).min(DUP_MAX + 1)
                    }
                }
                _ => -2,
            };
        }
        num
    }

    /// The character a token carries (glibc's `opr.c`).
    fn tok_char(&self) -> Option<u8> {
        match self.tok {
            Tok::Char(c) => Some(c),
            _ if self.tok_len == 2 => self.raw.get(self.at + 1).copied(),
            _ => self.byte(self.at),
        }
    }

    fn dup(&mut self, elem: Node) -> Result<Node, &'static str> {
        let (min, max) = match self.tok {
            Tok::OpenDup => {
                let mut start = self.number();
                if start == -1 {
                    if self.tok_char() == Some(b',') && matches!(self.tok, Tok::Char(_)) {
                        start = 0;
                    } else {
                        return Err("REG_BADBR");
                    }
                }
                let mut end = 0;
                if start != -2 {
                    end = if self.tok == Tok::CloseDup {
                        start
                    } else if self.tok_char() == Some(b',') && matches!(self.tok, Tok::Char(_)) {
                        self.number()
                    } else {
                        -2
                    };
                }
                if start == -2 || end == -2 {
                    return Err(if self.tok == Tok::End {
                        "REG_EBRACE"
                    } else {
                        "REG_BADBR"
                    });
                }
                if (end != -1 && start > end) || self.tok != Tok::CloseDup {
                    return Err("REG_BADBR");
                }
                if (if end == -1 { start } else { end }) > DUP_MAX {
                    return Err("REG_ESIZE");
                }
                (start as u32, (end >= 0).then_some(end as u32))
            }
            Tok::DupAsterisk => (0, None),
            Tok::DupPlus => (1, None),
            _ => (0, Some(1)),
        };
        self.fetch(false);
        if max == Some(0) {
            return Ok(Node::Empty);
        }
        Ok(Node::Repeat(Box::new(elem), min, max))
    }

    /// A bracket expression after its `[`.
    fn bracket(&mut self) -> Result<Node, &'static str> {
        self.at += self.tok_len;
        self.tok_len = 0;
        let mut set = [false; 256];
        let mut negate = false;
        let (mut t, mut len) = self.peek_bracket();
        if t == BTok::NonMatch {
            negate = true;
            self.at += len;
            (t, len) = self.peek_bracket();
        }
        if t == BTok::Close {
            t = BTok::Char(b']');
        }
        let mut first = true;
        loop {
            let start = self.bracket_element(t, len, first)?;
            first = false;
            (t, len) = self.peek_bracket();
            let mut range_end = None;
            if !matches!(start, Elem::Class(_) | Elem::Equiv(_)) {
                if t == BTok::End {
                    return Err("REG_EBRACK");
                }
                if t == BTok::Range {
                    self.at += len;
                    let (t2, len2) = self.peek_bracket();
                    if t2 == BTok::End {
                        return Err("REG_EBRACK");
                    }
                    if t2 == BTok::Close {
                        // A `-` before the closing bracket is ordinary.
                        self.at -= len;
                        t = BTok::Char(b'-');
                    } else {
                        range_end = Some((t2, len2));
                    }
                }
            }
            if let Some((t2, len2)) = range_end {
                let end = self.bracket_element(t2, len2, true)?;
                (t, len) = self.peek_bracket();
                let lo = range_point(&start)?;
                let hi = range_point(&end)?;
                if lo > hi {
                    return Err("REG_ERANGE");
                }
                for c in lo..=hi {
                    set[c as usize] = true;
                }
            } else {
                match start {
                    Elem::Byte(c) => set[c as usize] = true,
                    Elem::Coll(name) | Elem::Equiv(name) => {
                        if name.len() != 1 {
                            return Err("REG_ECOLLATE");
                        }
                        set[name[0] as usize] = true;
                    }
                    Elem::Class(name) => {
                        let f = class(&name, self.icase)?;
                        for c in 0..=255u8 {
                            if f(c) {
                                set[c as usize] = true;
                            }
                        }
                    }
                }
            }
            if t == BTok::End {
                return Err("REG_EBRACK");
            }
            if t == BTok::Close {
                break;
            }
        }
        self.at += len;
        if negate {
            for v in set.iter_mut() {
                *v = !*v;
            }
        }
        // The next `fetch` reads on from after the `]`.
        self.tok_len = 0;
        Ok(Node::Set(Box::new(set)))
    }

    fn peek_bracket(&self) -> (BTok, usize) {
        let Some(c) = self.byte(self.at) else {
            return (BTok::End, 0);
        };
        if c == b'[' {
            let c2 = self.byte(self.at + 1).unwrap_or(0);
            return match c2 {
                b'.' => (BTok::OpenColl, 2),
                b'=' => (BTok::OpenEquiv, 2),
                b':' => (BTok::OpenClass, 2),
                _ => (BTok::Char(c), 1),
            };
        }
        (
            match c {
                b'-' => BTok::Range,
                b']' => BTok::Close,
                b'^' => BTok::NonMatch,
                _ => BTok::Char(c),
            },
            1,
        )
    }

    fn bracket_element(
        &mut self,
        t: BTok,
        len: usize,
        accept_hyphen: bool,
    ) -> Result<Elem, &'static str> {
        self.at += len;
        match t {
            BTok::OpenColl | BTok::OpenEquiv | BTok::OpenClass => {
                let delim = match t {
                    BTok::OpenColl => b'.',
                    BTok::OpenEquiv => b'=',
                    _ => b':',
                };
                if self.at >= self.buf.len() {
                    return Err("REG_EBRACK");
                }
                let mut name = Vec::new();
                loop {
                    if name.len() >= 32 {
                        return Err("REG_EBRACK");
                    }
                    let ch = if t == BTok::OpenClass {
                        self.raw[self.at]
                    } else {
                        self.buf[self.at]
                    };
                    self.at += 1;
                    if self.at >= self.buf.len() {
                        return Err("REG_EBRACK");
                    }
                    if ch == delim && self.buf[self.at] == b']' {
                        break;
                    }
                    name.push(ch);
                }
                self.at += 1;
                Ok(match t {
                    BTok::OpenColl => Elem::Coll(name),
                    BTok::OpenEquiv => Elem::Equiv(name),
                    _ => Elem::Class(name),
                })
            }
            BTok::Range if !accept_hyphen => {
                // A `-` that does not start a range must end the list.
                if self.peek_bracket().0 != BTok::Close {
                    return Err("REG_ERANGE");
                }
                Ok(Elem::Byte(b'-'))
            }
            BTok::Range => Ok(Elem::Byte(b'-')),
            BTok::Close => Ok(Elem::Byte(b']')),
            BTok::NonMatch => Ok(Elem::Byte(b'^')),
            BTok::Char(c) => Ok(Elem::Byte(c)),
            BTok::End => Err("REG_EBRACK"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BTok {
    Char(u8),
    Range,
    Close,
    NonMatch,
    OpenColl,
    OpenEquiv,
    OpenClass,
    End,
}

enum Elem {
    Byte(u8),
    Coll(Vec<u8>),
    Equiv(Vec<u8>),
    Class(Vec<u8>),
}

fn range_point(e: &Elem) -> Result<u8, &'static str> {
    match e {
        Elem::Byte(c) => Ok(*c),
        Elem::Coll(name) if name.len() == 1 => Ok(name[0]),
        Elem::Coll(_) => Err("REG_ECOLLATE"),
        _ => Err("REG_ERANGE"),
    }
}

fn class_set(f: impl Fn(u8) -> bool, negate: bool) -> [bool; 256] {
    let mut set = [false; 256];
    for c in 0..=255u8 {
        set[c as usize] = f(c) != negate;
    }
    set
}

/// A `[:name:]` class in the C locale.
fn class(name: &[u8], icase: bool) -> Result<fn(u8) -> bool, &'static str> {
    let name = if icase && (name == b"upper" || name == b"lower") {
        b"alpha".as_slice()
    } else {
        name
    };
    Ok(match name {
        b"alnum" => |c: u8| c.is_ascii_alphanumeric(),
        b"cntrl" => |c: u8| c.is_ascii_control(),
        b"lower" => |c: u8| c.is_ascii_lowercase(),
        b"space" => |c: u8| crate::cnum::is_space(c),
        b"alpha" => |c: u8| c.is_ascii_alphabetic(),
        b"digit" => |c: u8| c.is_ascii_digit(),
        b"print" => |c: u8| (0x20..0x7f).contains(&c),
        b"upper" => |c: u8| c.is_ascii_uppercase(),
        b"blank" => |c: u8| c == b' ' || c == b'\t',
        b"graph" => |c: u8| (0x21..0x7f).contains(&c),
        b"punct" => |c: u8| c.is_ascii_punctuation(),
        b"xdigit" => |c: u8| c.is_ascii_hexdigit(),
        _ => return Err("REG_ECTYPE"),
    })
}

type Caps = Vec<Option<(usize, usize)>>;

fn has_backref(node: &Node) -> bool {
    match node {
        Node::BackRef(_) => true,
        Node::Group(n, _) | Node::Repeat(n, _, _) => has_backref(n),
        Node::Concat(v) | Node::Alt(v) => v.iter().any(has_backref),
        _ => false,
    }
}

/// With no branch or repetition, captures have one mandatory definition.
/// At most 32 AST nodes need fewer than 3*32 matcher/sequence frames, below
/// MATCH_DEPTH. Larger layouts retain the native fallback, not an error.
fn short_linear(node: &Node, count: &mut usize, captured: &mut u16) -> bool {
    if *count >= 32 {
        return false;
    }
    *count += 1;
    match node {
        Node::Alt(_) | Node::Repeat(_, _, _) => false,
        Node::Group(inner, g) => {
            if !short_linear(inner, count, captured) {
                return false;
            }
            if *g < 9 {
                *captured |= 1u16 << *g;
            }
            true
        }
        Node::BackRef(g) => *captured & (1u16 << *g) != 0,
        Node::Concat(nodes) => nodes.iter().all(|n| short_linear(n, count, captured)),
        _ => true,
    }
}

/// A bounded size estimate, not a bound on glibc's executed instructions.
fn expanded_size(node: &Node) -> usize {
    let cap = native::LIMIT + 1;
    match node {
        Node::Group(n, _) => expanded_size(n).saturating_add(2).min(cap),
        Node::Concat(v) | Node::Alt(v) => v.iter().fold(1usize, |n, child| {
            n.saturating_add(expanded_size(child)).min(cap)
        }),
        Node::Repeat(n, min, max) => {
            let copies = max.map_or(*min as usize + 1, |max| max as usize);
            expanded_size(n)
                .saturating_add(1)
                .saturating_mul(copies)
                .saturating_add(1)
                .min(cap)
        }
        _ => 1,
    }
}

impl Regex {
    /// `regcomp(pattern, REG_NOSUB | (icase ? REG_ICASE : 0))`; the error
    /// is glibc's error code name.
    pub(crate) fn new(pattern: &[u8], icase: bool) -> Result<Regex, &'static str> {
        let pattern = &pattern[..pattern
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(pattern.len())];
        let buf = if icase {
            pattern.to_ascii_uppercase()
        } else {
            pattern.to_vec()
        };
        let mut p = Parser {
            buf,
            raw: pattern,
            at: 0,
            tok: Tok::End,
            tok_len: 0,
            groups: 0,
            completed: 0,
            icase,
        };
        p.fetch(true);
        let root = p.reg_exp(0)?;
        let literal = Literal::new(&root);
        let native = if has_backref(&root) {
            let mut count = 0;
            let mut captured = 0;
            if short_linear(&root, &mut count, &mut captured) {
                None
            } else {
                Some(NativeRegex::new(pattern, icase, expanded_size(&root))?)
            }
        } else {
            None
        };
        Ok(Regex {
            root,
            icase,
            // BRE back references are single digits (\1 through \9).
            // Other groups affect structure but need no saved capture slots.
            groups: p.groups.min(9),
            literal,
            native,
        })
    }

    fn ch(&self, c: u8) -> u8 {
        if self.icase {
            c.to_ascii_uppercase()
        } else {
            c
        }
    }

    /// `regexec(regex, subject, 0, NULL, 0) == 0`.
    pub(crate) fn is_match(&self, s: &[u8]) -> bool {
        if let Some(native) = &self.native {
            return native.is_match(s);
        }
        if let Some(literal) = &self.literal {
            return literal.is_match(self, s);
        }
        BUDGET.with(|state| {
            state.set(Budget {
                work: MATCH_WORK,
                depth: 0,
                exhausted: false,
            })
        });
        let mut caps: Caps = vec![None; self.groups];
        let first_only = match &self.root {
            Node::Assert(Anchor::LineFirst | Anchor::BufFirst) => true,
            Node::Concat(nodes) => matches!(
                nodes.first(),
                Some(Node::Assert(Anchor::LineFirst | Anchor::BufFirst))
            ),
            _ => false,
        };
        // These exact root shapes fail their first assertion at every
        // positive start. Group/alternation structure is not inferred.
        let last_start = if first_only { 0 } else { s.len() };
        for start in 0..=last_start {
            if self.m(&self.root, s, start, &mut caps, &mut |_, _| true) {
                return true;
            }
            if exhausted() {
                break;
            }
        }
        false
    }

    fn assert(&self, a: Anchor, s: &[u8], i: usize) -> bool {
        let prev_word = i > 0 && is_word(s[i - 1]);
        let next_word = i < s.len() && is_word(s[i]);
        match a {
            Anchor::LineFirst | Anchor::BufFirst => i == 0,
            Anchor::LineLast | Anchor::BufLast => i == s.len(),
            Anchor::WordFirst => !prev_word && next_word,
            Anchor::WordLast => prev_word && !next_word,
            Anchor::WordDelim => prev_word != next_word,
            Anchor::NotWordDelim => prev_word == next_word,
        }
    }

    fn m(
        &self,
        n: &Node,
        s: &[u8],
        i: usize,
        caps: &mut Caps,
        k: &mut dyn FnMut(usize, &mut Caps) -> bool,
    ) -> bool {
        let Some(_frame) = Frame::enter() else {
            return false;
        };
        match n {
            Node::Empty => k(i, caps),
            Node::Char(c) => i < s.len() && self.ch(s[i]) == *c && k(i + 1, caps),
            Node::Any => i < s.len() && s[i] != 0 && k(i + 1, caps),
            Node::Set(set) => i < s.len() && set[self.ch(s[i]) as usize] && k(i + 1, caps),
            Node::Assert(a) => self.assert(*a, s, i) && k(i, caps),
            Node::Group(inner, g) => {
                let g = *g;
                if g >= self.groups {
                    return self.m(inner, s, i, caps, k);
                }
                self.m(inner, s, i, caps, &mut |j, caps: &mut Caps| {
                    let old = caps[g];
                    caps[g] = Some((i, j));
                    let r = k(j, caps);
                    caps[g] = old;
                    r
                })
            }
            Node::BackRef(g) => {
                let Some((a, b)) = caps[*g] else {
                    return false;
                };
                let len = b - a;
                i + len <= s.len()
                    && charge(len)
                    && (0..len).all(|d| self.ch(s[a + d]) == self.ch(s[i + d]))
                    && k(i + len, caps)
            }
            Node::Concat(v) => self.seq(v, 0, s, i, caps, k),
            Node::Alt(v) => v.iter().any(|alt| self.m(alt, s, i, caps, k)),
            Node::Repeat(inner, min, max) => self.rep(inner, *min, *max, 0, s, i, caps, k),
        }
    }

    fn seq(
        &self,
        v: &[Node],
        idx: usize,
        s: &[u8],
        i: usize,
        caps: &mut Caps,
        k: &mut dyn FnMut(usize, &mut Caps) -> bool,
    ) -> bool {
        let Some(_frame) = Frame::enter() else {
            return false;
        };
        if idx == v.len() {
            return k(i, caps);
        }
        self.m(&v[idx], s, i, caps, &mut |j, caps: &mut Caps| {
            self.seq(v, idx + 1, s, j, caps, k)
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn rep(
        &self,
        inner: &Node,
        min: u32,
        max: Option<u32>,
        count: u32,
        s: &[u8],
        i: usize,
        caps: &mut Caps,
        k: &mut dyn FnMut(usize, &mut Caps) -> bool,
    ) -> bool {
        let Some(_frame) = Frame::enter() else {
            return false;
        };
        let atomic = |c: u8| match inner {
            Node::Char(p) => self.ch(c) == *p,
            Node::Any => c != 0,
            Node::Set(set) => set[self.ch(c) as usize],
            _ => false,
        };
        if matches!(inner, Node::Char(_) | Node::Any | Node::Set(_)) {
            let limit = max.map_or(s.len() - i, |n| n as usize).min(s.len() - i);
            let mut n = 0;
            while n < limit && atomic(s[i + n]) {
                if !charge(1) {
                    return false;
                }
                n += 1;
            }
            if n < min as usize {
                return false;
            }
            for used in (min as usize..=n).rev() {
                if !charge(1) {
                    return false;
                }
                if k(i + used, caps) {
                    return true;
                }
                if exhausted() {
                    return false;
                }
            }
            return false;
        }
        // Zero-width iterations can update different capture slots. Continue
        // only while that state changes; all slots can change to this empty
        // position at most once along one branch, so this cannot cycle.
        let mut before = [None; 9];
        before[..caps.len()].copy_from_slice(caps);
        if max.is_none_or(|m| count < m)
            && self.m(inner, s, i, caps, &mut |j, caps: &mut Caps| {
                // Preserve new empty captures, then stop at their fixed point.
                if j == i && count >= min {
                    if caps.as_slice() != &before[..caps.len()] {
                        return self.rep(inner, min, max, count + 1, s, j, caps, k);
                    }
                    return k(i, caps);
                }
                self.rep(inner, min, max, count + 1, s, j, caps, k)
            })
        {
            return true;
        }
        count >= min && k(i, caps)
    }
}
