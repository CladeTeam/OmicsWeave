use super::corpus::*;
use crate::support::rng::Rng;

pub(super) fn header(rng: &mut Rng, names: &mut Vec<String>) -> Vec<Vec<u8>> {
    let mut lines: Vec<Vec<u8>> = Vec::new();
    let n_ref = rng.range(0, 4) as usize;
    let mut alts = ALT_NAMES.to_vec();
    for _ in 0..n_ref {
        let name = loop {
            let n = rng.pick(NAMES).to_string();
            if !names.contains(&n) {
                break n;
            }
        };
        let mut tags = vec![format!("SN:{name}"), format!("LN:{}", rng.pick(LENGTHS))];
        if rng.chance(20) {
            tags.swap(0, 1);
        }
        if rng.chance(20) {
            tags.push("M5:6aef897c3d6ff0c78aff06ac189178dd".into());
        }
        if rng.chance(20) && !alts.is_empty() {
            let k = rng.range(0, alts.len() as u64 - 1) as usize;
            let alt = alts.remove(k);
            tags.push(format!("AN:{alt},,{name}"));
            names.push(alt.to_string());
        }
        names.push(name);
        lines.push(format!("@SQ\t{}", tags.join("\t")).into_bytes());
    }
    let others: [(u64, &[u8]); 6] = [
        (40, b"@HD\tVN:1.6\tSO:coordinate"),
        (10, b"@HD\tVN:1.4"),
        (30, b"@RG\tID:rg1\tSM:sample one\tPL:ILLUMINA"),
        (
            30,
            b"@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem ref.fa r1.fq",
        ),
        (20, b"@PG\tID:samtools\tPN:samtools\tPP:bwa"),
        (30, b"@CO\tfree text\twith a tab"),
    ];
    for (percent, line) in others {
        if rng.chance(percent) {
            let at = rng.range(0, lines.len() as u64) as usize;
            lines.insert(at, line.to_vec());
        }
    }
    if !lines.is_empty() && rng.chance(10) {
        let at = rng.range(1, lines.len() as u64) as usize;
        lines.insert(at, b"@CO".to_vec());
    }
    if !lines.is_empty() && rng.chance(5) {
        let bad: &[&[u8]] = &[
            b"@SQ\tSN:nolength",
            b"@SQ\tSN:zero\tLN:0",
            b"@SQ\tSN:neg\tLN:-5",
            b"@SQ\tLN:10",
            b"@SQ\tSN:two\tLN:10\tLN:11",
            b"@SQ\tSN:same\tLN:10\tLN:10",
            b"@SQ\tSN:chr1\tLN:5\t",
            b"@SQ\tSNchr1\tLN:5",
            b"@RG\tSM:x",
            b"@PG\tPN:x",
            b"@XY\tAB:c",
            b"@HD",
            b"@Co\tx",
            b"@COx",
            b"@",
            b"@1Q\tSN:x",
        ];
        let at = rng.range(1, lines.len() as u64) as usize;
        lines.insert(at, rng.pick(bad).to_vec());
        if rng.chance(30) && !names.is_empty() {
            let dup = format!("@SQ\tSN:{}\tLN:9", names[0]);
            lines.push(dup.into_bytes());
        }
    }
    lines
}
