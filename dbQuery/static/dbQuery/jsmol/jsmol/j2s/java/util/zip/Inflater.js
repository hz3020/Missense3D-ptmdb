Clazz.declarePackage("java.util.zip");
Clazz.load(["com.jcraft.jzlib.Inflater"], "java.util.zip.Inflater", null, function(){
var c$ = Clazz.declareType(java.util.zip, "Inflater", com.jcraft.jzlib.Inflater);
Clazz.defineMethod(c$, "initialize", 
function(nowrap){
return this.init(0, nowrap);
}, "~B");
});
;//5.0.1-v2 Mon Jan 22 13:18:57 CST 2024
