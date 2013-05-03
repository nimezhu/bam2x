/**
 * Created with JetBrains WebStorm.
 * User: zhuxp
 * Date: 13-4-30
 * Time: 下午1:04
 * To change this template use File | Settings | File Templates.
 */

Array.prototype.max = function() {
    var max = this[0];
    var len = this.length;
    for (var i = 1; i < len; i++) if (this[i] > max) max = this[i];
    return max;
}
Array.prototype.min = function() {
    var min = this[0];
    var len = this.length;
    for (var i = 1; i < len; i++) if (this[i] < min) min = this[i];
    return min;
}


$(function() {
    Backbone.sync = function(method, model, success, error){
        success();
    }

    RNA = Backbone.Model.extend(
        { defaults: { id: "RNA" , struct2d:"(((...)))",seq:'' }}
    );

    RNAScore: Backbone.Model.extend(

    );


    RNAList = Backbone.Collection.extend({
        model: RNA,

        initialize: function(){

            this.bind("add", function( model ){

                var R2D = new RNA2DView({model:model});
                R2D.render();
            })



        }

    });



    LinksView = Backbone.View.extend(
    {   attributes:{},
        initialize: function(options) {
            if(options.attributes)
            {
                attributes=options.attributes;
            }
            _.bindAll(this,'render');
        },
        render: function()
        {
            var fill = d3.scale.ordinal()
                .domain(d3.range(4))
                .range(["#000000", "#FFDD89", "#957244", "#F26223"]);
            var stem_vector=parseStems(this.model.get("struct2d").trim());
            var struct2d=this.model.get("struct2d").trim();

            var len=this.model.get("struct2d").trim().length;
            var  link=this.el.select("#link");
            link.attr("transform","translate(250,250)");
            link.selectAll("g").remove();

            var link_stems=link.selectAll("g").data(stem_vector);
            var length=this.len;
            var link_g=link_stems.enter().append("g");

            link_g.append("path").attr("data",function(d) {console.log(d);return d})
                .attr("d",
                    d3.svg.chord()
                        .source(function(d,i) {console.log(d);return {startAngle:(d[0]-1)/len*6.20,
                            endAngle:d[1]/len*6.20}})
                        .target(function(d,i) {console.log(d); return {startAngle:(d[2]-1)/len*6.20,
                            endAngle:d[3]/len*6.20}})
                        .radius(150)


                )   .style("fill",  function(d,i){return fill(i%4)})

                .style("opacity", 0.5)
                .on("mouseover",function(d,i){$("#OUT").html(highlight_text(d)); d3.select(this).style("opacity",1.0);})
                .on("mouseout",function(d,i){$("#OUT").html(struct2d);d3.select(this).style("opacity",0.5)});

            function highlight_text(d)
            {
                var s=struct2d.split('');
                var state=0;
                var retv='';
                for (var i0=0;i0< s.length;i0++)
                {
                    if (i0+1==d[0]) {retv+='<span style="color: indianred">';}
                    if (i0+1==d[2]) {retv+='<span style="color: indianred">';}
                    retv+=s[i0];
                    if (i0+1==d[1]) {retv+='</span>';}
                    if (i0+1==d[3]) {retv+='</span>';}
                }
                return retv+d;
            }
        }
    }
    ),
    IdeogramView = Backbone.View.extend(
    {
       attributes: {},

       initialize: function(options) {
           if(options.attributes)
           {
               attributes=options.attributes;
           }
           _.bindAll(this,'render');

       },

       render: function()
       {
           //var ideogram=this.el.select("#ideogram");

           var ideogram=this.el.select("#ideogram");
           ideogram.selectAll("path").remove();
           ideogram.attr("transform","translate(250,250)");
           ideogram.append("path").attr("d", d3.svg.arc().outerRadius(220).innerRadius(210).startAngle(0).endAngle(6.20));
       }
    }
    )
    HighlightView = Backbone.View.extend(
    {
        attributes: {},

        initialize: function(options) {

            _.bindAll(this,'render');
            if(options.attributes)
            {
                attributes=options.attributes;
            }
            console.log("MODEL"+this.model);
            console.log("attributes"+this.attributes);
        },

        render: function()
        {
            console.log("attributes"+this.attributes);
            var stem_vector=parseStems(this.model.get("struct2d"));

            var length=this.model.get("struct2d").length;

            var highlight=this.el.select("#highlight");

            highlight.attr("transform","translate(250,250)");
            highlight.selectAll("g").remove();
            var stems=highlight.selectAll("g").data(stem_vector);
            var g_stems=stems.enter().append("g");

            g_stems.append("path").attr("fill","red")
                .attr("d",
                    d3.svg.arc()

                        .outerRadius(function(d) { return 160})
                        .innerRadius(function(d) { return 155})
                        .startAngle(function(d) { return (d[0]-1)/length*6.20})
                        .endAngle(function(d) { return d[1]/length*6.20})


                )
            g_stems.append("path").attr("fill","green")
                .attr("d",
                    d3.svg.arc()
                        .outerRadius(function(d) { return 160})
                        .innerRadius(function(d) { return 155})
                        .startAngle(function(d) { return (d[2]-1)/length*6.20})
                        .endAngle(function(d) { return d[3]/length*6.20})
                )

        }
    })


    PlotView = Backbone.View.extend(
        {


            attributes: {},
            initialize: function(options) {
                if(options.attributes)
                {
                    attributes=options.attributes;
                }

                _.bindAll(this,'render');

            },

            render: function()
            {
                var data=this.model.data;
                var len=this.model.data.length;
                console.log(data);
                var max_val=data.max();
                var min_val=data.min();
                console.log("min_val "+min_val);
                var min_radius=170; // 170 to 290 ?
                var max_radius=210;
                var plots=d3.select("#FIGURE").select("svg").select("#plot");
                plots.attr("transform","translate(250,250)");
                plots.selectAll("path").remove();
                var bars=plots.selectAll("path").data(data).enter().append("path");
                bars.attr("fill","#99FF44").attr("d",
                    d3.svg.arc()
                        .outerRadius(function(d) { return (d-min_val)/(max_val-min_val)*(max_radius-min_radius)+min_radius;})
                        .innerRadius(function(d) { return min_radius;}  )
                        .startAngle(function(d,i) { return i/len*6.20;})
                        .endAngle(function(d,i) {return (i+1)/len*6.20;})

                )


            }
        }



    );


    RNA2DView = Backbone.View.extend(
        {
        attributes: {"width":500,"height":500},
        initialize: function(options) {
            if(options.attributes)
            {
                attributes=options.attributes;
            }

        _.bindAll(this,'render');

         },

            render: function(  )
            {
                d3.select("#FIGURE").select("h2").text("FIGURE "+this.model.get("id"));

                var svg=d3.select("#FIGURE").select("svg").attr("width", this.attributes.width)   // <-- Here
                    .attr("height", this.attributes.height);


                var ideogramView = new IdeogramView({el:svg});

                var highlightView = new HighlightView({el:svg, model:this.model, attributes:this.attributes});
                var linksView = new LinksView({el:svg,model:this.model,attributes: this.attributes }) ;


                //var plotView = new PlotView({el:svg,model:test});
                ideogramView.render();
                highlightView.render();
                linksView.render();
                //plotView.render();




            }



        }
    )


    RNAView = Backbone.View.extend( {
      tagname:'li',
      events:{
        'click span.view' : 'view',
          'click span.delete': 'remove'
      },
            initialize: function() {

                _.bindAll(this,'render','view','unrender','remove');
                this.model.bind('change', this.render);
                this.model.bind('remove', this.unrender);

            },
            render: function()
            {
                $(this.el).html('<span style="color:black;">'+this.model.get("id")
                    +'</span> &nbsp; &nbsp; <span class="view" style="font-family:sans-serif; color:blue; cursor:pointer;">[view]</span> '
                    +'<span class="delete" style="cursor:pointer; color:red; font-family:sans-serif;">[delete]</span>');
                return this;
            },
            view: function()
            {
                var R2D = new RNA2DView({model:this.model});
                R2D.render();

                $("#OUT").text(this.model.get("struct2d"));//
                return this;
            },
            unrender: function(){
                $(this.el).remove();

            },
            remove: function(){
                this.model.destroy();

            }

        }





    );


    RNAListView = Backbone.View.extend(
       {   el: $('#OUTPUTLIST'), // el attaches to existing element
           attributes: {"width":500,"height":500},
           events: {
               'click button#add': 'addRNA',
               'click button#reset': 'reset',
               'click button#clear' : 'clear',
               'click button#plotScore' : 'plotScore'
           },


          initialize: function()
          {
              _.bindAll(this, 'render', 'addRNA', 'appendRNA');
             this.collection= new RNAList();
               // every function that uses 'this' as the current object should be in here


              this.collection.bind('add', this.appendRNA); // collection event binder

              this.counter = 0;
              this.render();
           },
           render: function()
           {

           },
           appendRNA: function(item)
           {
               var itemView = new RNAView({
                   model: item
               });

               $('ul', this.el).append(itemView.render().el);

           },
           addRNA: function(){
               this.counter++;

               var local_ss = $('#INPUT').val().trim();


               var name="RNA_"+this.counter.toString();

               var item = new RNA({id:name,struct2d:local_ss});
               //var itemScore = new RNAScore({RNA:item,score:{"data":local_val}});
               console.log(name+" "+local_ss);
               console.log(this.collection)
               this.collection.add(item);

           },

           reset: function() {
               var count = this.collection.size();
               for (var i = count-1; i > -1; i--) {
                   this.collection.remove(this.collection.at(i));
               }
               $("#OUTPUTLIST").text('');
           },
           clear: function()
           {
             $("#INPUT").val("");

           },



           plotScore: function()
           {
               var data=$("#PLOT").val().trim().split(",");
              // var svg=d3.select("#FIGURE").select("svg").attr("width", this.attributes.width)   // <-- Here
              //     .attr("height", this.attributes.height);
              // var p=new PlotView({model:{'el':svg,'data':data}});
               var p=new PlotView({model:{'data':data}});
               p.render();
           }
       }
    )







    var view = new RNAListView({el: 'body'});
    //var scoreView = new PlotControl();
});