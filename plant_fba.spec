/*
A KBase module: plant_fba
*/

module plant_fba {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;
    
    /*
    @optional input_columns
    */
    typedef structure {
	string input_ws;
	string input_expression_matrix;
	string input_fbamodel;
	string input_columns;
	string output_reaction_matrix;
    } IntegrateAbundancesParams;

    funcdef integrate_abundances_with_metabolism(IntegrateAbundancesParams input_params)
        returns (ReportResults output_report)
	authentication required;

    typedef structure {
	string input_ws;
	string input_genome;

	string output_ws;
	string output_fbamodel;

	string template;
	string template_ws;
    } ReconstructMetabolismParams;

    funcdef reconstruct_plant_metabolism(ReconstructMetabolismParams input_params)
        returns (ReportResults output_report)
	authentication required;
};
