from bokeh.layouts import row, grid
from bokeh.models import Slope, Select, ColumnDataSource, CustomJS, NumeralTickFormatter
from bokeh.plotting import figure, output_file, save

import time
import json
import math

from plant_fba.core.fetch_plantseed_impl import FetchPlantSEEDImpl
#from fetch_plantseed_impl import FetchPlantSEEDImpl

class GenerateFigureImpl:

    def log(self,message, prefix_newline=False):
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
        print(('\n' if prefix_newline else '') + time_str + ': ' + message)

    def __init__(self):
        pass

    def generate_figure(self, columns, category_select=None, 
                          genome_features=None, reaction_scores=None, 
                          reaction_percentiles=None):

        # To include with figure object
        TOOLTIPS = [
            ("reaction", "@tooltip"),
            ("(x,y)", "($x, $y)")
            ]

        ##################################################################
        # The output figure will be saved using the 'grid' function
        # Each row in the figure will be from a pair of columns in the matrix
        # The first scatterplot will be the general "genome-features background"
        # The second scatterplot will be the reaction percentiles and subsystems
        # The third "column" in a row will contain the subsystem select

        figure_grid = list()
        for first_column in range(len(columns)):
            for second_column in range(len(columns)):
                if(first_column <= second_column):
                    continue

                # Row of figures for the pair of conditions
                figure_row = list()

                ##################################################################
                # For the first scatterplot, it is optional
                if( genome_features is not None or reaction_scores is not None ):

                    # Find range for axes
                    x_max = math.ceil(max(genome_features[columns[first_column]]))
                    y_max = math.ceil(max(genome_features[columns[second_column]]))
                    plot_max = max([x_max,y_max])

                    bokeh_fig = figure(x_range = (0.0,plot_max),
                                       y_range = (0.0,plot_max))
                    bokeh_fig.xaxis.axis_label = columns[first_column]
                    bokeh_fig.yaxis.axis_label = columns[second_column]
                    bokeh_fig.title.text = "Genome Features Expression Abundances"

                    genome_source = ColumnDataSource(data=dict(genome_features))

                    # Plot as black and visible
                    scatter_fig = bokeh_fig.circle(x=columns[0], y=columns[1], source=genome_source,
                                                   color='black', size=4, visible=True)
                    
                    reaction_source = ColumnDataSource(data=dict(reaction_scores))

                    # Plot as red
                    scatter_fig = bokeh_fig.circle(x=columns[0], y=columns[1], source=reaction_source,
                                                   color='red', size=6, visible=True)

                    slope_line = Slope(gradient=1, y_intercept=0, line_color="red")
                    bokeh_fig.add_layout(slope_line)

                    figure_row.append(bokeh_fig)

                ##################################################################
                # For the second scatterplot
                if( reaction_percentiles is not None ):

                    ##################################################################
                    # Set up parent figure object

                    bokeh_fig = figure(tooltips=TOOLTIPS,
                                       x_range = (0.0,1.0),
                                       y_range = (0.0,1.0))
                    bokeh_fig.xaxis.axis_label = columns[first_column]
                    bokeh_fig.yaxis.axis_label = columns[second_column]
                    bokeh_fig.xaxis.formatter = NumeralTickFormatter(format="0.0")
                    bokeh_fig.yaxis.formatter = NumeralTickFormatter(format="0.0")
                    bokeh_fig.title.text = "Model Reactions Percentile Rank (p<0.01)"

                    ##################################################################
                    # The data is transformed into ColumnDataSource object to allow for CustomJS to work
                    # The source_dict stores the data after it's been transformed into ColumnDataSource
                    # The scatter_dict stores the individual bokeh scatterplots for rendering in CustomJS
                    source_dict=dict()
                    scatter_dict=dict()

                    ##################################################################
                    # For the background data, all the data is captured under a single 'All' key
                    # It is added first, so that it will always be in the background
                    # It is intentionally made visible and won't be changed in the CustomJS
                    # Transform
                    source = ColumnDataSource(data=dict(reaction_percentiles['All']))
                    # Store transformation
                    source_dict['All']=source
                    # Plot as black and visible
                    scatter_fig = bokeh_fig.circle(x=columns[0], y=columns[1], source=source,
                                                   color='color', size='size', fill_alpha='fill_alpha',
                                                   visible=True)
                    # Store plot
                    scatter_dict['All']=scatter_fig

                    ##################################################################
                    # For the foreground data, the scatter plot for each subsystem is create
                    # separately, but made invisible, to be used with the Select dropdown

                    for scatter in reaction_percentiles.keys():
                        # Not using the 'All' background data
                        if(scatter == 'All'):
                            continue

                        # Transform
                        source = ColumnDataSource(data=dict(reaction_percentiles[scatter]))
                        # Store transformation
                        source_dict[scatter]=source
                        # Plot as red but not visible
                        scatter_fig = bokeh_fig.circle(x=columns[0], y=columns[1], source=source,
                                                       color='color', size='size', fill_alpha='fill_alpha',
                                                       visible=False)
                        # Store plot
                        scatter_dict[scatter]=scatter_fig

                    # Add red central slope
                    slope_line = Slope(gradient=1, y_intercept=0, line_color="red")
                    bokeh_fig.add_layout(slope_line)

                    # Add parent figure to row of figures
                    figure_row.append(bokeh_fig)

                    # Add subsystem selector
                    # Starts with default value of "None" and allows user to pick one
                    # whereupon, according to CustomJS code below, it'll become visible
                    subsystem_select = Select(title="Select Subsystem:", 
                                              value="None", options=['None']+sorted(category_select))

                    # Add JS callback
                    callback = CustomJS(args=dict(source=source_dict,
                                                  figs=scatter_dict,
                                                  subsystem_select=subsystem_select),
                                        code="""
console.log("Updating")
for (let scatter in source){

    // Only choose subsystem
    if(scatter == 'All'){
        continue
    }

    // Chosen subsystem
    if(scatter == subsystem_select.value){
        
        figs[scatter].visible=true

        // Iterate through datapoints to make sure they are red and of a larger size
        for (let i = 0; i < source[scatter].data['color'].length; i++) {
            // This is where I would scale with p-value
            source[scatter].data['color'][i] = 'red' 
            source[scatter].data['size'][i] = 8
        }

    } else {

        // Here we have to make sure that the non-chosen subsystems are not visible
        figs[scatter].visible=false

        // default values, but this is really un-necessary
        for (let i = 0; i < source[scatter].data['color'].length; i++) {
            source[scatter].data['color'][i] = 'black'
            source[scatter].data['size'][i] = 6
        }
    }
    // Actually show change in plot
    source[scatter].change.emit()
}
""")
                    subsystem_select.js_on_change('value', callback)

                    # Add subsystem selector to row of figures
                    figure_row.append(subsystem_select)

                    # Add row of figures to grid
                    figure_grid.append(figure_row)

        return figure_grid

def main():

    plantseed = FetchPlantSEEDImpl()
    reactions_data = plantseed.fetch_reactions()
    
    # Fetch ReactionMatrix
    with open("../../..//data/Ath_H13_Reaction_Matrix.json") as fh:
        reaction_matrix = json.load(fh)

    conditions = reaction_matrix['data']['col_ids']
    reactions = reaction_matrix['data']['row_ids']

    # The data is re-formatted so that it can be translated into a ColumnDataSource (CDS) object
    # So each reaction has a row, wherein it has the x and y values (headed by the conditions
    # they came from) and values for 'color', 'size', and a 'tooltip', which allows us to use
    # the CDS object to change them in the CustomJS callback

    reaction_scores_dict=dict()
    subsystem_select_list=list()
    for i in range(len(reactions)):

        ##################################################################
        # Set-up reaction-specific data for each datum in the plot
        # Generate subsystem attribute for tooltip
        rxn = reactions[i].split('_')[0]
        ss_string = ",".join(reactions_data[rxn]['subsystems'])
        for ss in reactions_data[rxn]['subsystems']:

            # Set up list of subsystems for selection
            if(ss not in subsystem_select_list):
                subsystem_select_list.append(ss)

            # if(ss not in reaction_scores_dict):
            #    reaction_scores_dict[ss]=dict()
            #    reaction_scores_dict[ss]['color']=list()
            #    reaction_scores_dict[ss]['size']=list()
            #    reaction_scores_dict[ss]['tooltip']=list()

            # default values for plot
            # reaction_scores_dict[ss]['color'].append('red')
            # reaction_scores_dict[ss]['size'].append(8)
            # reaction_scores_dict[ss]['tooltip'].append(reactions[i]+" ("+ss_string+")")

        if('All' not in reaction_scores_dict):
            reaction_scores_dict['All']=dict()
            reaction_scores_dict['All']['color']=list()
            reaction_scores_dict['All']['size']=list()
            reaction_scores_dict['All']['tooltip']=list()

        # default values for plot
        reaction_scores_dict['All']['color'].append('black')
        reaction_scores_dict['All']['size'].append(6)
        reaction_scores_dict['All']['tooltip'].append(reactions[i]+" ("+ss_string+")")

        ##################################################################
        # Iterate conditions and save the x,y value for each reaction

        for j in range(len(conditions)):

            # Some reactions are literally not mapped because of a lack of data (rather than zero)
            # but I force these to zero to try and keep things consistent in the plot which
            # auto-scales
            str_value = "{0:.2f}".format(reaction_matrix['data']['values'][i][j])
            if(float(str_value) < -1.00):
                str_value = "0.00"

            # save the values under 'All' for background data
            if(conditions[j] not in reaction_scores_dict['All']):
                reaction_scores_dict['All'][conditions[j]]=list()

            reaction_scores_dict['All'][conditions[j]].append(float(str_value))

            # save the values for each subsystem
            # for ss in reactions_data[rxn]['subsystems']:
            #    if(conditions[j] not in reaction_scores_dict[ss]):
            #        reaction_scores_dict[ss][conditions[j]]=list()
            #    reaction_scores_dict[ss][conditions[j]].append(float(str_value))

    loaded_scores_dict = dict()
    with open("shit.json") as fh:
        loaded_scores_dict = json.load(fh)

    subsystem_select_list = ["None"]
    figure_generator = GenerateFigureImpl()
    figure_grid = figure_generator.generate_figure(conditions, category_select=subsystem_select_list,
                                                   reaction_percentiles=loaded_scores_dict)

    ##################################################################
    # Save output
    file_name='test_multiscatter_output.html'
    output_file(file_name)
    save(grid(figure_grid))

if(__name__ == "__main__"):
    main()
