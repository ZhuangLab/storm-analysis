import storm_analysis


def test_frc():
	
	in_list = storm_analysis.get_data("test/data/test_drift_mlist.bin")
	results = storm_analysis.get_path_output_test("test_drift_frc.txt")

	from storm_analysis.frc.frc_calc2d import calc2d
	calc2d(in_list, results, False)


if __name__ == "__main__":
	test_frc()