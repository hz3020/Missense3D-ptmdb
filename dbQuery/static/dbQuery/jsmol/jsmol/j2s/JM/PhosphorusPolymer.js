Clazz.declarePackage("JM");
Clazz.load(["JM.BioPolymer"], "JM.PhosphorusPolymer", null, function(){
var c$ = Clazz.declareType(JM, "PhosphorusPolymer", JM.BioPolymer);
Clazz.makeConstructor(c$, 
function(monomers){
Clazz.superConstructor (this, JM.PhosphorusPolymer, []);
this.set(monomers);
this.hasStructure = true;
}, "~A");
});
;//5.0.1-v2 Sun Feb 04 10:11:05 CST 2024
