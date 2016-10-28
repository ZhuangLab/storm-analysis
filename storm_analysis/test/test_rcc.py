import storm_analysis


def test_rcc():
	
	bin = storm_analysis.get_data("test/data/test_drift_mlist.bin")
	drift = storm_analysis.get_path_output_test("test_drift.txt")

	from storm_analysis.rcc.rcc_drift_correction import rcc
	rcc(bin, drift, 2000, 1)


if __name__ == "__main__":
	test_rcc()