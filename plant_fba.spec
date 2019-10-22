/*
A KBase module: plant_fba
*/

module plant_fba {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;
    
    typedef structure {
	string input_ws;
	string input_expression_matrix;
	string input_fbamodel;
	string output_reaction_matrix;
    } IntegrateAbundancesParams;

    funcdef integrate_abundances_with_metabolism(IntegrateAbundancesParams input_params)
        returns (ReportResults output_report)
	authentication required;
};
