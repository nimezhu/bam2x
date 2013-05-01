/**
 * Created with JetBrains WebStorm.
 * User: zhuxp
 * Date: 13-4-30
 * Time: 下午1:04
 * To change this template use File | Settings | File Templates.
 */




$(function() {
    Backbone.sync = function(method, model, success, error){
        success();
    }

    RNA = Backbone.Model.extend(
        { defaults: { id: "RNA" , struct2d:"(((...)))",seq:""}
    }

    );

    RNAList = Backbone.Collection.extend({
        model: RNA,

        initialize: function(){

            this.bind("add", function( model ){

                var R2D = new RNA2DView(model);
                R2D.render();
            })



        }

    });

    RNA2DView = Backbone.View.extend(
        {
        initialize: function(model) {
            this.model=model;
        //_.bindAll(this,'render')

    },
    render: function()
    {
        this.render_figure_header(this.model.get("id"));
        this.render_figure(this.model.get("struct2d"));
        return this;
    },
      //need revise to delete the #part


            render_figure_header: function( id )
            {
                d3.select("#FIGURE").select("h2").text("FIGURE "+id);
            }
            ,
            render_figure: function( struct2d )
            {

                var fill = d3.scale.ordinal()
                    .domain(d3.range(4))
                    .range(["#000000", "#FFDD89", "#957244", "#F26223"]);
                var stem_vector=parseStems(struct2d);
                console.log("stem len:"+stem_vector.length);
                console.log("struct string:" + struct2d);
                var w=500;
                var h=500;
                var svg=d3.select("#FIGURE").select("svg").attr("width", w)   // <-- Here
                    .attr("height", h);
                var ideogram=svg.select("#ideogram");
                ideogram.attr("transform","translate(250,250)");
                ideogram.append("path").attr("d", d3.svg.arc().outerRadius(220).innerRadius(210).startAngle(0).endAngle(6.20));

                var highlight=svg.select("#highlight")
                highlight.attr("transform","translate(250,250)");
                highlight.selectAll("g").remove();
                var stems=highlight.selectAll("g").data(stem_vector);

                var length=struct2d.length;

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


                var  link=svg.select("#link");
                link.attr("transform","translate(250,250)");
                link.selectAll("g").remove();

                var link_stems=link.selectAll("g").data(stem_vector);
                var link_g=link_stems.enter().append("g");

                link_g.append("path").attr("data",function(d) {return d})
                    .attr("d",
                        d3.svg.chord()
                            .source(function(d,i) {return {startAngle:(stem_vector[i][0]-1)/length*6.20,
                                endAngle:stem_vector[i][1]/length*6.20}})
                            .target(function(d,i) {return {startAngle:(stem_vector[i][2]-1)/length*6.20,
                                endAngle:stem_vector[i][3]/length*6.20}})
                            .radius(150)


                    )   .style("fill",  function(d,i){return fill(i%4)})

                    .style("opacity", 0.5)
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
                var R2D = new RNA2DView(this.model);
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
       {   el: $('OUTPUTLIST'), // el attaches to existing element
           events: {
               'click button#add': 'addRNA',
               'click button#reset': 'reset',
               'click button#clear' : 'clear'
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


               //$(this.el).append("<ul></ul>");
              /*
               _(this.collection.models).each(function(item){ // in case collection is not empty
                   self.appendRNA(item);
               }, this);
                */

               console.log("LEN:"+this.collection.length);
           },
           appendRNA: function(item)
           {
               var itemView = new RNAView({
                   model: item
               });

               $('ul', this.el).append(itemView.render().el);
               //$().append(itemView.render().el);
           },
           addRNA: function(){
               this.counter++;

               var local_ss = $('#INPUT').val();

               var name="RNA_"+this.counter.toString();

               var item = new RNA({id:name,struct2d:local_ss});
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

           }

       }
    )





    var view = new RNAListView({el: 'body'});
});