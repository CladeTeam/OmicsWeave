use omicsweave_convert::Budget;

fn budget() -> Budget {
    Budget::new(2, 32 << 20)
}

mod metadata;
mod options;
mod planning;
mod publication;
