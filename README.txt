#########################################################
#
# author: Yoni Leibner
# description: Analyzer class allows to easy visualization
#              and exploration of neurons from multiple
#              point of view, and also to create neuron
#              movies.
# date of modification: 16.11.2022
#
#########################################################

instalation:
    clone the rep from github
    open console
    go the folder directory
    install requirements.txt: pip install -r Neuron_analysis_tool/requirements.txt
    install Neuron_analysis_tool: pip install Neuron_analysis_tool

explanation:

    this project helps the analysis on neuron to equivelant cables/ dendorgams and voltage decays
    the package have 2 main part:

    peeler - a gui on matplotlib to help peel expinantinal decays from short pulse (see Rall ref), a word document will be attached shortly
    analyzer - this tool can help you anlyze your neuron!!!! it have:
                morphology plot - where you can color code each segment
                morphology value plot - where you can give a value for color coding the morphology
                dendogram plot - where you can color code each segment in both physical units and electrical units
                dendogram value plot - where you can give a value for color coding the morphology
                attenuation plot - ploting the attenuation along the tree -> y axis can be normed, and x axis (distance from origin) can be scaled to electrical units
                video creator - create a video on voltage/channel conductances of preset protocol or protocol of your choosing


please look on the 3 jupyter notebook in this repo and for any question contact me on yoni.leibner@mail.huji.ac.il


