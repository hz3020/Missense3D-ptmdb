Clazz.declarePackage("J.adapter.readers.quantum");
Clazz.load(["J.adapter.smarter.AtomSetCollectionReader"], "J.adapter.readers.quantum.GaussianWfnReader", null, function(){
var c$ = Clazz.declareType(J.adapter.readers.quantum, "GaussianWfnReader", J.adapter.smarter.AtomSetCollectionReader);
Clazz.overrideMethod(c$, "initializeReader", 
function(){
this.continuing = false;
});
});
;//5.0.1-v2 Sun Feb 04 10:11:05 CST 2024
