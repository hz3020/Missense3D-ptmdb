Clazz.declarePackage("JS");
Clazz.load(["java.util.Hashtable", "JU.P3", "$.V3"], "JS.WyckoffFinder", ["JU.JSJSONParser", "$.M4", "$.Measure", "$.P4", "$.PT", "$.Rdr", "$.SB", "JV.FileManager"], function(){
var c$ = Clazz.decorateAsClass(function(){
this.positions = null;
this.npos = 0;
this.ncent = 0;
this.centerings = null;
this.centeringStr = null;
Clazz.instantialize(this, arguments);}, JS, "WyckoffFinder", null);
/*LV!1824 unnec constructor*/Clazz.makeConstructor(c$, 
function(map){
if (map != null) {
var wpos = map.get("wpos");
this.positions = wpos.get("pos");
this.npos = this.positions.size();
var cent = wpos.get("cent");
if (cent != null) {
this.ncent = cent.size();
this.centeringStr =  new Array(this.ncent);
this.centerings =  new Array(this.ncent);
for (var i = this.ncent; --i >= 0; ) {
var s = cent.get(i);
this.centeringStr[i] = s;
this.centerings[i] = JS.WyckoffFinder.toPoint(s);
}
}}}, "java.util.Map");
Clazz.defineMethod(c$, "getWyckoffFinder", 
function(vwr, sgname){
var helper = JS.WyckoffFinder.helpers.get(sgname);
if (helper == null) {
helper = JS.WyckoffFinder.createHelper(this, vwr, sgname);
}if (helper == null) {
if (JS.WyckoffFinder.nullHelper == null) JS.WyckoffFinder.nullHelper =  new JS.WyckoffFinder(null);
JS.WyckoffFinder.helpers.put(sgname, JS.WyckoffFinder.nullHelper);
}return helper;
}, "JV.Viewer,~S");
Clazz.defineMethod(c$, "findPositionFor", 
function(p, letter){
if (this.positions != null) {
for (var i = this.npos; --i >= 0; ) {
var map = this.positions.get(i);
if (map.get("label").equals(letter)) {
var coords = map.get("coord");
if (coords != null) JS.WyckoffFinder.getWyckoffCoord(coords, 0, letter).project(p);
return p;
}}
}return null;
}, "JU.P3,~S");
Clazz.defineMethod(c$, "getStringInfo", 
function(uc, p, returnType){
var info = this.createInfo(uc, p, returnType);
return (info == null ? "?" : info);
}, "JS.UnitCell,JU.P3,~N");
c$.toPoint = Clazz.defineMethod(c$, "toPoint", 
function(xyz){
var s = JU.PT.split(xyz, ",");
return JU.P3.new3(JU.PT.parseFloatFraction(s[0]), JU.PT.parseFloatFraction(s[1]), JU.PT.parseFloatFraction(s[2]));
}, "~S");
c$.wrap = Clazz.defineMethod(c$, "wrap", 
function(xyz, sb){
return sb.appendC('(').append(xyz).appendC(')');
}, "~S,JU.SB");
Clazz.defineMethod(c$, "createInfo", 
function(uc, p, returnType){
switch (returnType) {
case 42:
var sb =  new JU.SB();
this.getCenteringStr(-1, sb);
for (var i = this.npos; --i >= 0; ) {
var map = this.positions.get(i);
var label = map.get("label");
sb.appendC('\n').append(label);
if (i == 0) {
sb.append(" (x,y,z)");
} else {
JS.WyckoffFinder.getList(map.get("coord"), label, sb);
}}
return sb.toString();
case -1:
case -2:
case -3:
for (var i = this.npos; --i >= 0; ) {
var map = this.positions.get(i);
var label = map.get("label");
if (i == 0) {
switch (returnType) {
case -1:
return label;
case -2:
return "(x,y,z)";
case -3:
return map.get("label") + "  (x,y,z)";
}
}var coords = map.get("coord");
for (var c = 0, n = coords.size(); c < n; c++) {
var coord = JS.WyckoffFinder.getWyckoffCoord(coords, c, label);
if (coord.contains(this, uc, p)) {
switch (returnType) {
case -1:
return label;
case -2:
return coord.asString(null, true).toString();
case -3:
var sbc =  new JU.SB();
sbc.append(label).appendC(' ');
this.getCenteringStr(-1, sbc).appendC(' ');
JS.WyckoffFinder.getList(coords, label, sbc);
return sbc.toString();
}
}}
}
break;
default:
var letter = "" + String.fromCharCode(returnType);
for (var i = this.npos; --i >= 0; ) {
var map = this.positions.get(i);
if (map.get("label").equals(letter)) {
return (i == 0 ? "(x,y,z)" : JS.WyckoffFinder.getList(map.get("coord"), letter, null).toString());
}}
break;
}
return null;
}, "JS.UnitCell,JU.P3,~N");
c$.createHelper = Clazz.defineMethod(c$, "createHelper", 
function(w, vwr, sgname){
var itno = JU.PT.parseInt(JU.PT.split(sgname, ":")[0]);
if (itno >= 1 && itno <= 230) {
var resource = JS.WyckoffFinder.getResource(w, vwr, "ita_" + itno + ".json");
if (resource != null) {
var its = resource.get("its");
if (its != null) {
for (var i = its.size(); --i >= 0; ) {
var map = its.get(i);
if (sgname.equals(map.get("itaFull"))) {
var helper =  new JS.WyckoffFinder(map);
JS.WyckoffFinder.helpers.put(sgname, helper);
return helper;
}}
}}}return null;
}, "~O,JV.Viewer,~S");
Clazz.defineMethod(c$, "getCenteringStr", 
function(index, sb){
if (sb == null) sb =  new JU.SB();
if (this.ncent == 0) return sb;
if (index >= 0) {
sb.appendC('+');
return JS.WyckoffFinder.wrap(this.centeringStr[index], sb);
}for (var i = 0; i < this.ncent; i++) {
sb.appendC('+');
JS.WyckoffFinder.wrap(this.centeringStr[i], sb);
}
return sb;
}, "~N,JU.SB");
c$.getList = Clazz.defineMethod(c$, "getList", 
function(coords, letter, sb){
if (sb == null) sb =  new JU.SB();
for (var c = 0, n = coords.size(); c < n; c++) {
var coord = JS.WyckoffFinder.getWyckoffCoord(coords, c, letter);
sb.append(" ");
coord.asString(sb, false);
}
return sb;
}, "JU.Lst,~S,JU.SB");
c$.getWyckoffCoord = Clazz.defineMethod(c$, "getWyckoffCoord", 
function(coords, c, label){
var coord = coords.get(c);
if ((typeof(coord)=='string')) {
coords.set(c, coord =  new JS.WyckoffFinder.WyckoffCoord(coord, label));
}return coord;
}, "JU.Lst,~N,~S");
c$.getResource = Clazz.defineMethod(c$, "getResource", 
function(w, vwr, resource){
try {
var r = JV.FileManager.getBufferedReaderForResource(vwr, w, "JS/", "sg/json/" + resource);
var data =  new Array(1);
if (JU.Rdr.readAllAsString(r, 2147483647, false, data, 0)) {
return  new JU.JSJSONParser().parse(data[0], true);
}} catch (e) {
System.err.println(e.getMessage());
}
return null;
}, "~O,JV.Viewer,~S");
/*if3*/;(function(){
var c$ = Clazz.decorateAsClass(function(){
this.type = 0;
this.xyz = null;
this.label = null;
this.thisCentering = "";
this.op = null;
this.point = null;
this.line = null;
this.plane = null;
Clazz.instantialize(this, arguments);}, JS.WyckoffFinder, "WyckoffCoord", null);
Clazz.makeConstructor(c$, 
function(xyz, label){
this.xyz = xyz;
this.label = label;
this.create(xyz);
}, "~S,~S");
Clazz.defineMethod(c$, "contains", 
function(w, uc, p){
var slop = uc.getPrecision();
this.thisCentering = null;
if (this.checkLatticePt(p, slop)) return true;
if (w.centerings == null) return false;
for (var i = w.centerings.length; --i >= 0; ) {
JS.WyckoffFinder.WyckoffCoord.pc.add2(p, w.centerings[i]);
uc.unitize(JS.WyckoffFinder.WyckoffCoord.pc);
if (this.checkLatticePt(JS.WyckoffFinder.WyckoffCoord.pc, slop)) {
this.thisCentering = w.centeringStr[i];
return true;
}}
return false;
}, "JS.WyckoffFinder,JS.UnitCell,JU.P3");
Clazz.defineMethod(c$, "project", 
function(p){
switch (this.type) {
case 1:
p.setT(this.point);
break;
case 2:
JU.Measure.projectOntoAxis(p, this.point, this.line, JS.WyckoffFinder.WyckoffCoord.vt);
break;
case 3:
JU.Measure.getPlaneProjection(p, this.plane, JS.WyckoffFinder.WyckoffCoord.vt, JS.WyckoffFinder.WyckoffCoord.vt);
p.setT(JS.WyckoffFinder.WyckoffCoord.vt);
break;
}
}, "JU.P3");
Clazz.defineMethod(c$, "asString", 
function(sb, withCentering){
if (sb == null) sb =  new JU.SB();
JS.WyckoffFinder.wrap(this.xyz, sb);
if (withCentering && this.thisCentering != null) {
sb.appendC('+');
JS.WyckoffFinder.wrap(this.thisCentering, sb);
}return sb;
}, "JU.SB,~B");
Clazz.defineMethod(c$, "checkLatticePt", 
function(p, slop){
if (this.checkPoint(p, slop)) return true;
for (var z = 62, i = -2; i < 3; i++) {
for (var j = -2; j < 3; j++) {
for (var k = -2; k < 3; k++, z--) {
if (z == 0) continue;
JS.WyckoffFinder.WyckoffCoord.p3.set(i, j, k);
JS.WyckoffFinder.WyckoffCoord.p3.add(p);
if (this.checkPoint(JS.WyckoffFinder.WyckoffCoord.p3, slop)) {
System.out.println(this.label + " " + this.xyz + " found for " + i + " " + j + " " + k);
return true;
}}
}
}
return false;
}, "JU.P3,~N");
Clazz.defineMethod(c$, "checkPoint", 
function(p, slop){
var d = 1;
switch (this.type) {
case 1:
d = this.point.distance(p);
break;
case 2:
JS.WyckoffFinder.WyckoffCoord.p1.setT(p);
JU.Measure.projectOntoAxis(JS.WyckoffFinder.WyckoffCoord.p1, this.point, this.line, JS.WyckoffFinder.WyckoffCoord.vt);
d = JS.WyckoffFinder.WyckoffCoord.p1.distance(p);
break;
case 3:
d = Math.abs(JU.Measure.getPlaneProjection(p, this.plane, JS.WyckoffFinder.WyckoffCoord.vt, JS.WyckoffFinder.WyckoffCoord.vt));
break;
}
if (d < slop) System.out.println("success! " + this.label + " " + this.xyz + " " + d + " " + p);
return d < slop;
}, "JU.P3,~N");
Clazz.defineMethod(c$, "create", 
function(p){
var nxyz = (p.indexOf('x') >= 0 ? 1 : 0) + (p.indexOf('y') >= 0 ? 1 : 0) + (p.indexOf('z') >= 0 ? 1 : 0);
if (nxyz == 1 || nxyz == 2) {
var a =  Clazz.newFloatArray (16, 0);
var v = JU.PT.split(this.xyz, ",");
JS.WyckoffFinder.WyckoffCoord.getRow(v[0], a, 0);
JS.WyckoffFinder.WyckoffCoord.getRow(v[1], a, 4);
JS.WyckoffFinder.WyckoffCoord.getRow(v[2], a, 8);
a[15] = 1;
this.op = JU.M4.newA16(a);
}switch (nxyz) {
case 0:
this.type = 1;
this.point = JS.WyckoffFinder.toPoint(p);
break;
case 1:
this.type = 2;
JS.WyckoffFinder.WyckoffCoord.p1.set(0.19, 0.53, 0.71);
this.op.rotTrans(JS.WyckoffFinder.WyckoffCoord.p1);
JS.WyckoffFinder.WyckoffCoord.p2.set(0.51, 0.27, 0.64);
this.op.rotTrans(JS.WyckoffFinder.WyckoffCoord.p2);
JS.WyckoffFinder.WyckoffCoord.p2.sub2(JS.WyckoffFinder.WyckoffCoord.p2, JS.WyckoffFinder.WyckoffCoord.p1);
JS.WyckoffFinder.WyckoffCoord.p2.normalize();
this.point = JU.P3.newP(JS.WyckoffFinder.WyckoffCoord.p1);
this.line = JU.V3.newV(JS.WyckoffFinder.WyckoffCoord.p2);
break;
case 2:
this.type = 3;
JS.WyckoffFinder.WyckoffCoord.p1.set(0.19, 0.51, 0.73);
this.op.rotTrans(JS.WyckoffFinder.WyckoffCoord.p1);
JS.WyckoffFinder.WyckoffCoord.p2.set(0.23, 0.47, 0.86);
this.op.rotTrans(JS.WyckoffFinder.WyckoffCoord.p2);
JS.WyckoffFinder.WyckoffCoord.p3.set(0.1, 0.2, 0.);
this.op.rotTrans(JS.WyckoffFinder.WyckoffCoord.p3);
this.plane = JU.Measure.getPlaneThroughPoints(JS.WyckoffFinder.WyckoffCoord.p1, JS.WyckoffFinder.WyckoffCoord.p2, JS.WyckoffFinder.WyckoffCoord.p3, null, null,  new JU.P4());
break;
case 3:
break;
}
}, "~S");
c$.getRow = Clazz.defineMethod(c$, "getRow", 
function(s, a, rowpt){
s = JU.PT.rep(s, "-", "+-");
s = JU.PT.rep(s, "x", "*x");
s = JU.PT.rep(s, "y", "*y");
s = JU.PT.rep(s, "z", "*z");
s = JU.PT.rep(s, "-*", "-");
s = JU.PT.rep(s, "+*", "+");
var part = JU.PT.split(s, "+");
for (var p = part.length; --p >= 0; ) {
s = part[p];
if (s.length == 0) continue;
var pt = 3;
if (s.indexOf('.') >= 0) {
var d = JU.PT.parseFloat(s);
a[rowpt + pt] = d;
continue;
}var i0 = 0;
var sgn = 1;
switch ((s.charAt(0)).charCodeAt(0)) {
case 45:
sgn = -1;
case 42:
i0++;
break;
}
var v = 0;
for (var i = s.length, f2 = 0; --i >= i0; ) {
var c = s.charAt(i);
switch ((c).charCodeAt(0)) {
case 120:
pt = 0;
v = 1;
break;
case 121:
pt = 1;
v = 1;
break;
case 122:
pt = 2;
v = 1;
break;
case 47:
f2 = 1;
v = 1 / v;
case 42:
sgn *= v;
v = 0;
break;
default:
var u = "0123456789".indexOf(c);
if (u < 0) System.err.println("WH ????");
if (v == 0) {
v = u;
} else {
f2 = (f2 == 0 ? 10 : f2 * 10);
v += f2 * u;
}break;
}
}
a[rowpt + pt] = sgn * v;
}
}, "~S,~A,~N");
Clazz.overrideMethod(c$, "toString", 
function(){
return this.asString(null, false).toString();
});
c$.p1 =  new JU.P3();
c$.p2 =  new JU.P3();
c$.p3 =  new JU.P3();
c$.pc =  new JU.P3();
c$.vt =  new JU.V3();
/*eoif3*/})();
c$.nullHelper = null;
c$.helpers =  new java.util.Hashtable();
});
;//5.0.1-v2 Mon Feb 05 08:36:38 CST 2024
