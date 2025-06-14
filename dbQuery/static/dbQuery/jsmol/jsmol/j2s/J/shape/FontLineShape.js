Clazz.declarePackage("J.shape");
Clazz.load(["J.shape.Shape"], "J.shape.FontLineShape", null, function(){
var c$ = Clazz.decorateAsClass(function(){
this.tickInfos = null;
this.font3d = null;
Clazz.instantialize(this, arguments);}, J.shape, "FontLineShape", J.shape.Shape);
Clazz.prepareFields (c$, function(){
this.tickInfos =  new Array(4);
});
Clazz.overrideMethod(c$, "initShape", 
function(){
this.translucentAllowed = false;
});
Clazz.defineMethod(c$, "setPropFLS", 
function(propertyName, value){
if ("tickInfo" === propertyName) {
var t = value;
if (t.ticks == null) {
if (t.type.equals(" ")) {
this.tickInfos[0] = this.tickInfos[1] = this.tickInfos[2] = this.tickInfos[3] = null;
return;
}t = null;
}this.tickInfos["xyz".indexOf(t.type) + 1] = t;
return;
}if ("font" === propertyName) {
this.font3d = value;
return;
}}, "~S,~O");
Clazz.overrideMethod(c$, "getShapeState", 
function(){
return null;
});
});
;//5.0.1-v2 Sun Feb 04 10:11:05 CST 2024
