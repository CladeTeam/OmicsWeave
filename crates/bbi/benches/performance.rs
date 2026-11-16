mod performance {
    pub mod schema;
    pub mod variable_blocks;
    pub mod writer;
}
use performance::{schema, variable_blocks, writer};
#[path = "../tests/support/allocator.rs"]
mod allocator;
fn main() {
    writer::writer_benchmark();
    writer::warm_flush_allocations();
    variable_blocks::variable_block_benchmark();
    schema::auto_sql_benchmark();
}

#[path = "../tests/support/writer_fixture.rs"]
mod writer_fixture;
