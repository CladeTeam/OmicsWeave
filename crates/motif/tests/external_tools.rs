//! Optional comparisons against MEME 5.5.9. All tests are explicitly ignored.
//! Set each OMICSWEAVE_* variable and run with --ignored; missing tools fail.
#[allow(dead_code)]
#[path = "support/glam2.rs"]
mod glam2_support;
#[allow(dead_code)]
#[path = "support/mcast.rs"]
mod mcast_support;
#[path = "support/rng.rs"]
mod rng;

#[path = "external/fimo.rs"]
mod fimo;
#[path = "external/glam2.rs"]
mod glam2;
#[path = "external/mast.rs"]
mod mast;
#[path = "external/mcast.rs"]
mod mcast;
