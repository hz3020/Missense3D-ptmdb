Clazz.declarePackage("java.net");
(function(){
var c$ = Clazz.declareType(java.net, "URLEncoder", null);
c$.encode = Clazz.defineMethod(c$, "encode", 
function(s){
return encodeURIComponent(s);
}, "~S");
})();
;//5.0.1-v2 Mon Jan 22 13:18:57 CST 2024
