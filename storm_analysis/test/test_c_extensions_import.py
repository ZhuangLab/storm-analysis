# Test every extensions have been correctly build

import platform

def test_c_extensions_import():
	"Try to import all C extension"
	
	import storm_analysis.sa_library.ia_utilities_c
	import storm_analysis.sa_library.multi_fit_c
	import storm_analysis.sa_library.grid_c
	import storm_analysis.sa_library.matched_filter_c

	import storm_analysis.frc.frc_c

	import storm_analysis.fista.fista_decon_utilities_c
	import storm_analysis.fista.fista_fft_c

	import storm_analysis.fista.fista_fft_c

	import storm_analysis.dbscan.dbscan_c

	import storm_analysis.decon_storm.mlem_c

	import storm_analysis.sCMOS.scmos_utilities_c

	import storm_analysis.simulator.zernike_c
	import storm_analysis.simulator.drawgaussians

	import storm_analysis.spliner.cubic_spline_c
	import storm_analysis.spliner.cubic_fit_c

	import storm_analysis.rolling_ball_bgr.rolling_ball_lib_c

	import storm_analysis.sa_utilities.fitz_c
	import storm_analysis.sa_utilities.apply_drift_correction_c
	import storm_analysis.sa_utilities.avemlist_c
	import storm_analysis.sa_utilities.tracker_c

        #import storm_analysis.L1H.homotopy_imagea_c

	if platform.system() == 'Windows':
                import storm_analysis.L1H.fista_lib_c


if __name__ == "__main__":
	test_c_extensions_import()
