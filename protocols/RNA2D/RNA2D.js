/**
 * Created with JetBrains WebStorm.
 * User: zhuxp
 * Date: 13-4-21
 * Time: 上午9:33
 * To change this template use File | Settings | File Templates.
 */


function parseStems(ss)
{
    var sList=ss.split('');
    //console.log(ss);
    //console.log(sList);
    var level=0;
    var posRegister=Array();
    var pos=0;
    var pairRegister=Array();
    var pairs=Array();
    var stems=Array();
    for(var x in sList)
    {
        //console.log(x);
        //console.log(sList[x]);
        processChar(sList[x]);
        //console.log(level);

    }
    if(pairRegister.length>0)
    {
      dumpPairRegister();
    }

   // console.log(pairs[0])  ;
    for(var i in stems) console.log(stems[i]);

    return stems;

    function processChar(a)
    {
        pos+=1;
        if (a=="(")
        {
             posRegister.push(pos);

        }
        if (a==")")
        {   var first=posRegister.pop();
            var pair=[first,pos] ;
            pairs.push(pair);
            if(pairRegister.length>0)
            {
              var k=pairRegister.length-1;
              if((pairRegister[k][0]-first>1) || (pos-pairRegister[k][1] > 1))
              {
                  dumpPairRegister();

              }
            }
            pairRegister.push(pair);
            //console.log(pair);


        }
        if (a==".")
        {

        }

    }


    function dumpPairRegister()
    {
        var last=pairRegister.length-1;
        stems.push([pairRegister[last][0],pairRegister[0][0],pairRegister[0][1],pairRegister[last][1]]);
        pairRegister=[];
    }
}



