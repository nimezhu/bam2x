from __future__ import print_function
import os
import sys
import logging
import argparse
from  bam2x import IO,TableIO,DBI

template = '''
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>circos plot result</title>
  <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js"></script>
  <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>
  <script src="http://underscorejs.org/underscore-min.js"></script>
<script src="http://backbonejs.org/backbone-min.js"></script>
<script src="http://v.zhu.land/gist/circos/circos.js"></script>
  <style>
  body {
  font:10px sans-serif;
  }
  </style>
</head>

<body>
<div id="canvas">
</div>

</body>

<script> 
     
var colors = d3.scale.category10().domain(d3.range(0,10));



function plot(data) {
      d3.select('#canvas').text('');
      var outerRadius=250;
      var innteRadius=70;
      var plotHeight=30;
      var bedHeight=10;
      //var cx = outerRadius + 30;
      var cy = outerRadius + 30;
      //var cy = document.getElementById('canvas').clientHeight/2
      var cx = document.getElementById('canvas').clientWidth/2;
      var gapHeight= 5
      var nowRadius=outerRadius;
      var svg = d3.select("#canvas").append("svg").attr("id","svg");
      svg.attr("width",cx*2)
      svg.attr("height",document.getElementById('canvas').clientHeight)
      var collection = new IdeogramsCollection();
      var ideograms = []
      for (var i in data.ideograms){
          ideograms.push(new IdeogramModel(data.ideograms[i]))
      }
      collection.add(ideograms);
      var ideogramView = new IdeogramView({});
      var ideogramTrack = new IdeogramTrack({"collection":collection,"el":svg.append("g"),"cx":cx,"cy":cy,"outerRadius":nowRadius,"innerRadius":nowRadius-bedHeight});
      ideogramTrack.render(true);
     
     nowRadius=nowRadius-bedHeight-gapHeight;
     for( var i in data.tracks)
     {
     track=data.tracks[i];
     var plots=[];
     if (track.type=="plot")
     {
     for( var j in track.values)
     {
         var model=new PlotModel(track.values[j]);
         if (track.color)
         {
             model.color=track.color
         }
         plots.push(model);

     }
      var plotsCollection= new PlotsCollection(plots,{"name":track.name});
      var plotTrack = new PlotTrack({"collection":plotsCollection,"el":svg.append("g"),"cx":cx,"cy":cy,'outerRadius':nowRadius,'innerRadius':nowRadius-plotHeight});
      plotTrack.render(ideogramTrack);
      nowRadius-=plotHeight+gapHeight
    };

    if ( track.type=="bedgraph"){
        var collection= new BedGraphsCollection()
        collection.add(track.values)
        console.log(collection)
        var bedGraphTrack = new BedGraphView({"collection":collection,"el":svg.append("g"),"cx":cx,"cy":cy,"outerRadius":nowRadius,"innerRadius":nowRadius-plotHeight});
        bedGraphTrack.render(ideogramTrack);
        nowRadius-=plotHeight+gapHeight
        }
        

    if (track.type=="links")
    {
        var links = new LinksCollection();

        for(var i in track.values){
            links.push(new LinkModel(track.values[i]));
         }
        var linkTrack = new LinkTrack({"collection":links,"el":svg.append("g"),"cx":cx,"cy":cy,'radius':nowRadius});
        linkTrack.render(ideogramTrack);
    };



    }
}

var data = $DATA

$(document).ready(function() {
    plot(data);
    window.onresize = function(evt) {plot(data)}
});


</script>
</body>
'''
def help():
    return "json to circos html"
def set_parser(parser):
    #parser.add_argument("-m",type=str,choices=("seq","cDNA","cdna","cds","utr5","utr3"),dest="method")
    pass    
    
def run(args):
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    json=fin.read()
    a=template.replace("$DATA",json)
    print(a,file=out)

if __name__=="__main__":
    from bam2x.IO import parser_factory
    p=parser_factory(description=help())
    set_parser(p)
    run(p.parse_args())







