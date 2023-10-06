from simreaduntil.usecase_helpers.readfish_wrappers import NanoSimMapper, replace_ru_mapper


def test_alignment():
    list(NanoSimMapper("dummy").map_reads_2([(None, "chr11-NC-000011_76599_perfect_proc0:0_F_0_967_0", "A"*1000, 1000, None)]))
    
    with replace_ru_mapper(replace=True) as mapper_class:
        assert mapper_class is NanoSimMapper
    with replace_ru_mapper(replace=False) as mapper_class:
        assert mapper_class is not NanoSimMapper