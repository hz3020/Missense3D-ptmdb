Clazz.declarePackage("J.shapebio");
Clazz.load(["J.shapebio.Rockets"], "J.shapebio.Cartoon", null, function(){
var c$ = Clazz.declareType(J.shapebio, "Cartoon", J.shapebio.Rockets);
Clazz.defineMethod(c$, "initShape", 
function(){
Clazz.superCall(this, J.shapebio.Cartoon, "initShape", []);
this.madDnaRna = 1000;
});
});
;//5.0.1-v2 Sun Feb 04 10:11:05 CST 2024
