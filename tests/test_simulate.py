import pytest

import data

@pytest.mark.parametrize(
    "genome_fn,model,readlen,fraglen,fragsd",
    [
        (data.ref, "HS20", 100, 450, 150)
    ])
def test_RandomReadGenerator(genome_fn, model, readlen, fraglen, fragsd):
    from metacov.simulate import generate_reads
    reads = generate_reads(genome_fn, model, readlen, fraglen, fragsd)
    with reads as (genome_len, fqpair):
        for n in range(10):
            next(fqpair)
            assert fqpair.rlen == readlen+1
