Clazz.declarePackage("JU");
(function(){
var c$ = Clazz.decorateAsClass(function(){
this.x = 0;
this.y = 0;
this.width = 0;
this.height = 0;
Clazz.instantialize(this, arguments);}, JU, "Rectangle", null);
Clazz.defineMethod(c$, "contains", 
function(X, Y){
return (X >= this.x && Y >= this.y && X - this.x < this.width && Y - this.y < this.height);
}, "~N,~N");
})();
;//5.0.1-v2 Sun Feb 04 10:11:05 CST 2024
