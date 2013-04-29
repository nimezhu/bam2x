/**
 * Created with JetBrains WebStorm.
 * User: zhuxp
 * Date: 13-4-26
 * Time: 下午4:40
 * To change this template use File | Settings | File Templates.
 */



$(function() {



    RNAsList = Backbone.Collection.extend({
        initialize: function(){
            this.bind("add", function( model ){

                view.render( model );
            })
        }
    });

    RNAView = Backbone.View.extend({

        tagName: 'li',

        events: {
            'click #add-input':  'getRNA'
        },

        initialize: function() {
            this.RNAslist = new RNAsList;
            _.bindAll(this, 'render');
        },

        getRNA: function() {
            var local_ss = $('#RNA2D').val();
            this.RNAslist.add( {ss: local_ss} );
        },

        render: function( model ) {
            $("#OUTPUTLIST").append("<li>"+ model.get("ss"));
            a=parseStems(model.get("ss"));
            $("#OUTPUTLIST").append("<br>");
            for (var i in a)
            {
                $("#OUTPUTLIST").append(a[i]+"<br>");
            }
            $("#OUTPUTLIST").append("</li>");
            console.log('rendered')
        }

    });

    var view = new RNAView({el: 'body'});
});