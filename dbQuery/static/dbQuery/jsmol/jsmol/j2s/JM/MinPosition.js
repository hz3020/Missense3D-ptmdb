Clazz.declarePackage("JM");
Clazz.load(["JM.MinObject"], "JM.MinPosition", null, function(){
var c$ = Clazz.declareType(JM, "MinPosition", JM.MinObject);
Clazz.makeConstructor(c$, 
function(data, ddata){
Clazz.superConstructor (this, JM.MinPosition, []);
this.data = data;
this.ddata = ddata;
}, "~A,~A");
});
;//5.0.1-v2 Sun Feb 04 10:11:05 CST 2024
