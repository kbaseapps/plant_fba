# plant_fba release notes
=========================================

## v1.1.0
* Major update to the output of both Apps

    * "Reconstruct Plant Metabolism" has a new output table which
      captures the relationship between the functional annotation in
      the Genome and the reactions in the resulting FBAModel.

    * "Integrate Abundances with Metabolism" has a new scatterplot
      that complements the information in the table by allowing users
      to select and highlight the reaction expression scores for
      individual subsystems. We also expanded the table show the
      "enzyme-limiting" gene from which the reaction expression score
      is derived.

## v1.0.1
* Fixing of a bug where it assumes that the _functions_ field of a
genome feature object is present

## v1.0.0
* First major release of plant_fba
* Contains reconstruct_plant_metabolism function (exposed as
"Reconstruct Plant Metabolism" App)
* Contains integrate_abundances_with_metabolism (exposed as "Integrate
Abundances with Metabolism" App)

## v0.0.0
* Module created by kb-sdk init
