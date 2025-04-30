Clazz.declarePackage("JM");
Clazz.load(["JM.BioPolymer"], "JM.CarbohydratePolymer", null, function(){
var c$ = Clazz.declareType(JM, "CarbohydratePolymer", JM.BioPolymer);
Clazz.makeConstructor(c$, 
function(monomers){
Clazz.superConstructor (this, JM.CarbohydratePolymer, []);
this.set(monomers);
this.type = 3;
}, "~A");
});
;//5.0.1-v2 Sun Feb 04 10:11:05 CST 2024
