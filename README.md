heatmap-app
===========

A Shiny App for Interactive Gene Expression HeatMap Visualisation

This is a small Shiny App to display an interactive clustering/heatmap of differential gene expression. It takes as
data a list of DGELRT objects created from the edgeR package. Each slot of the list is name after the experiment
(the list names will show up at the bottom of the heatmap)

This is work in progress where hopefull I will be able to add table of genes from generated from the heatmap. The 
goal would be to click on a dendrogram node and have the associated genes from that branch populating the tables, one 
tab per experiment.

Fork me and help us develop that tools!
