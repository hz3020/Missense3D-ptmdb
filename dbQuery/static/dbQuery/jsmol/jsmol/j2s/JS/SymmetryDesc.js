Clazz.declarePackage("JS");
Clazz.load(["JU.P3", "$.V3"], "JS.SymmetryDesc", ["java.util.Hashtable", "JU.BS", "$.Lst", "$.M4", "$.Measure", "$.P4", "$.PT", "$.Quat", "$.SB", "JS.T", "JS.SpaceGroup", "$.Symmetry", "$.SymmetryOperation", "JU.Escape", "$.Logger"], function(){
var c$ = Clazz.decorateAsClass(function(){
this.modelSet = null;
Clazz.instantialize(this, arguments);}, JS, "SymmetryDesc", null);
/*LV!1824 unnec constructor*/Clazz.defineMethod(c$, "set", 
function(modelSet){
this.modelSet = modelSet;
return this;
}, "JM.ModelSet");
c$.getType = Clazz.defineMethod(c$, "getType", 
function(id){
var type;
if (id == null) return 1073742001;
if (id.equalsIgnoreCase("matrix")) return 12;
if (id.equalsIgnoreCase("description")) return 1825200146;
if (id.equalsIgnoreCase("axispoint")) return 134217751;
if (id.equalsIgnoreCase("time")) return 268441089;
if (id.equalsIgnoreCase("info")) return 1275068418;
if (id.equalsIgnoreCase("element")) return 1086326789;
if (id.equalsIgnoreCase("invariant")) return 36868;
type = JS.T.getTokFromName(id);
if (type != 0) return type;
type = JS.SymmetryDesc.getKeyType(id);
return (type < 0 ? type : 1073742327);
}, "~S");
c$.getKeyType = Clazz.defineMethod(c$, "getKeyType", 
function(id){
if ("type".equals(id)) id = "_type";
for (var type = 0; type < JS.SymmetryDesc.keys.length; type++) if (id.equalsIgnoreCase(JS.SymmetryDesc.keys[type])) return -1 - type;

return 0;
}, "~S");
c$.nullReturn = Clazz.defineMethod(c$, "nullReturn", 
function(type){
switch (type) {
case 135176:
return ";draw ID sym* delete;draw ID sg* delete;";
case 1073741961:
case 1825200146:
case 1073741974:
case 1145047050:
case 1145047053:
case 11:
case 1073742078:
return "";
case 1153433601:
return  new JU.BS();
default:
return null;
}
}, "~N");
c$.getInfo = Clazz.defineMethod(c$, "getInfo", 
function(io, type){
if (io.length == 0) return "";
if (type < 0 && -type <= JS.SymmetryDesc.keys.length && -type <= io.length) return io[-1 - type];
switch (type) {
case 1073742327:
case 1073741982:
return io;
case 1275068418:
var lst =  new java.util.Hashtable();
for (var j = 0, n = io.length; j < n; j++) {
var key = (j == 3 ? "draw" : j == 7 ? "axispoint" : JS.SymmetryDesc.keys[j]);
if (io[j] != null) lst.put(key, io[j]);
}
return lst;
case 1073741961:
return io[0] + "  \t" + io[2];
case 1145047050:
return io[0];
case 1145047053:
return io[19];
case 1073742078:
return io[1];
default:
case 1825200146:
return io[2];
case 134217764:
if (!JU.Logger.debugging) return io[3];
case 135176:
return io[3] + "\nprint " + JU.PT.esc(io[0] + " " + io[2]);
case 1145047051:
return io[4];
case 1073742178:
return io[5];
case 12289:
return io[6];
case 134217751:
return io[7];
case 1073741854:
return io[8];
case 134217729:
return io[9];
case 12:
return io[10];
case 1814695966:
return io[11];
case 4160:
return io[12];
case 268441089:
return io[13];
case 134217750:
return io[14];
case 1140850696:
return io[15];
case 1073741974:
return io[16];
case 1086326789:
return  Clazz.newArray(-1, [io[6], io[7], io[8], io[14], io[5]]);
case 36868:
return (io[6] != null ? io[6] : io[8] != null ?  Clazz.newArray(-1, [io[7], io[8], io[5]]) : io[5] != null ? "none" : io[14] != null ? io[14] : "identity");
}
}, "~A,~N");
c$.getInfoBS = Clazz.defineMethod(c$, "getInfoBS", 
function(type){
var bsInfo =  new JU.BS();
if (type < 0 && -type <= JS.SymmetryDesc.keys.length) {
bsInfo.set(-1 - type);
return bsInfo;
}switch (type) {
case 0:
case 1153433601:
case 1073742001:
case 1073742327:
case 1073741982:
case 1275068418:
bsInfo.setBits(0, JS.SymmetryDesc.keys.length);
break;
case 1073741961:
bsInfo.set(0);
bsInfo.set(2);
break;
case 1145047050:
bsInfo.set(0);
break;
case 1145047053:
bsInfo.set(19);
break;
case 1073742078:
bsInfo.set(1);
break;
default:
case 1825200146:
bsInfo.set(2);
break;
case 135176:
bsInfo.set(0);
bsInfo.set(2);
bsInfo.set(3);
break;
case 1145047051:
bsInfo.set(4);
break;
case 1073742178:
bsInfo.set(5);
break;
case 12289:
bsInfo.set(6);
break;
case 134217751:
bsInfo.set(7);
break;
case 1073741854:
bsInfo.set(8);
break;
case 134217729:
bsInfo.set(9);
break;
case 12:
bsInfo.set(10);
break;
case 1814695966:
bsInfo.set(11);
break;
case 4160:
bsInfo.set(12);
break;
case 268441089:
bsInfo.set(13);
break;
case 134217750:
bsInfo.set(14);
break;
case 1140850696:
bsInfo.set(15);
break;
case 1073741974:
bsInfo.set(16);
break;
case 1086326789:
case 36868:
bsInfo.set(5);
bsInfo.set(6);
bsInfo.set(7);
bsInfo.set(8);
bsInfo.set(14);
bsInfo.set(22);
break;
}
return bsInfo;
}, "~N");
Clazz.defineMethod(c$, "createInfoArray", 
function(op, uc, ptFrom, ptTarget, id, scaleFactor, options, haveTranslation, bsInfo, isSpaceGroup, isSpaceGroupAll){
if (!op.isFinalized) op.doFinalize();
var matrixOnly = (bsInfo.cardinality() == 1 && bsInfo.get(10));
var isTimeReversed = (op.timeReversal == -1);
if (scaleFactor == 0) scaleFactor = 1;
JS.SymmetryDesc.ptemp.set(0, 0, 0);
JS.SymmetryDesc.vtrans.set(0, 0, 0);
var plane = null;
var pta00 = (ptFrom == null || Float.isNaN(ptFrom.x) ?  new JU.P3() : ptFrom);
if (ptTarget != null) {
JS.SymmetryDesc.pta01.setT(pta00);
JS.SymmetryDesc.pta02.setT(ptTarget);
uc.toFractional(JS.SymmetryDesc.pta01, true);
uc.toFractional(JS.SymmetryDesc.pta02, true);
op.rotTrans(JS.SymmetryDesc.pta01);
JS.SymmetryDesc.ptemp.setT(JS.SymmetryDesc.pta01);
uc.unitize(JS.SymmetryDesc.pta01);
JS.SymmetryDesc.vtrans.setT(JS.SymmetryDesc.pta02);
uc.unitize(JS.SymmetryDesc.pta02);
if (JS.SymmetryDesc.pta01.distanceSquared(JS.SymmetryDesc.pta02) >= 1.96E-6) return null;
JS.SymmetryDesc.vtrans.sub(JS.SymmetryDesc.ptemp);
}var m2 = JU.M4.newM4(op);
m2.add(JS.SymmetryDesc.vtrans);
if (bsInfo.get(10) && ptTarget != null && pta00.equals(ptTarget)) {
m2.m00 = Math.round(m2.m00);
m2.m01 = Math.round(m2.m01);
m2.m02 = Math.round(m2.m02);
m2.m03 = Math.round(m2.m03);
m2.m10 = Math.round(m2.m10);
m2.m11 = Math.round(m2.m11);
m2.m12 = Math.round(m2.m12);
m2.m13 = Math.round(m2.m13);
m2.m20 = Math.round(m2.m20);
m2.m21 = Math.round(m2.m21);
m2.m22 = Math.round(m2.m22);
m2.m23 = Math.round(m2.m23);
}var isMagnetic = (op.timeReversal != 0);
if (matrixOnly && !isMagnetic) {
var im = JS.SymmetryDesc.getKeyType("matrix");
var o =  new Array(-im);
o[-1 - im] = m2;
return o;
}var ftrans =  new JU.V3();
JS.SymmetryDesc.pta01.set(1, 0, 0);
JS.SymmetryDesc.pta02.set(0, 1, 0);
var pta03 = JU.P3.new3(0, 0, 1);
JS.SymmetryDesc.pta01.add(pta00);
JS.SymmetryDesc.pta02.add(pta00);
pta03.add(pta00);
var pt0 = JS.SymmetryDesc.rotTransCart(op, uc, pta00, JS.SymmetryDesc.vtrans);
var pt1 = JS.SymmetryDesc.rotTransCart(op, uc, JS.SymmetryDesc.pta01, JS.SymmetryDesc.vtrans);
var pt2 = JS.SymmetryDesc.rotTransCart(op, uc, JS.SymmetryDesc.pta02, JS.SymmetryDesc.vtrans);
var pt3 = JS.SymmetryDesc.rotTransCart(op, uc, pta03, JS.SymmetryDesc.vtrans);
var vt1 = JU.V3.newVsub(pt1, pt0);
var vt2 = JU.V3.newVsub(pt2, pt0);
var vt3 = JU.V3.newVsub(pt3, pt0);
JS.SymmetryOperation.approx6Pt(JS.SymmetryDesc.vtrans);
JS.SymmetryDesc.vtemp.cross(vt1, vt2);
var haveInversion = (JS.SymmetryDesc.vtemp.dot(vt3) < 0);
if (haveInversion) {
pt1.sub2(pt0, vt1);
pt2.sub2(pt0, vt2);
pt3.sub2(pt0, vt3);
}var q = JU.Quat.getQuaternionFrame(pt0, pt1, pt2).div(JU.Quat.getQuaternionFrame(pta00, JS.SymmetryDesc.pta01, JS.SymmetryDesc.pta02));
var qF = JU.Quat.new4(q.q1, q.q2, q.q3, q.q0);
var info = JU.Measure.computeHelicalAxis(pta00, pt0, qF);
var pa1 = JU.P3.newP(info[0]);
var ax1 = JU.P3.newP(info[1]);
var ang1 = Clazz.floatToInt(Math.abs(JU.PT.approx((info[3]).x, 1)));
var pitch1 = JS.SymmetryOperation.approx((info[3]).y);
if (haveInversion) {
pt1.add2(pt0, vt1);
pt2.add2(pt0, vt2);
pt3.add2(pt0, vt3);
}var trans = JU.V3.newVsub(pt0, pta00);
if (trans.length() < 0.1) trans = null;
var ptinv = null;
var ipt = null;
var ptref = null;
var w = 0;
var margin = 0;
var isTranslation = (ang1 == 0);
var isRotation = !isTranslation;
var isInversionOnly = false;
var isMirrorPlane = false;
var isTranslationOnly = !isRotation && !haveInversion;
if (isRotation || haveInversion) {
trans = null;
}if (haveInversion && isTranslation) {
ipt = JU.P3.newP(pta00);
ipt.add(pt0);
ipt.scale(0.5);
ptinv = pt0;
isInversionOnly = true;
} else if (haveInversion) {
var d = (pitch1 == 0 ?  new JU.V3() : ax1);
var f = 0;
switch (ang1) {
case 60:
f = 0.6666667;
break;
case 120:
f = 2;
break;
case 90:
f = 1;
break;
case 180:
ptref = JU.P3.newP(pta00);
ptref.add(d);
pa1.scaleAdd2(0.5, d, pta00);
if (ptref.distance(pt0) > 0.1) {
trans = JU.V3.newVsub(pt0, ptref);
JS.SymmetryDesc.setFractional(uc, trans, JS.SymmetryDesc.ptemp, null);
ftrans.setT(JS.SymmetryDesc.ptemp);
} else {
trans = null;
}JS.SymmetryDesc.vtemp.setT(ax1);
JS.SymmetryDesc.vtemp.normalize();
w = -JS.SymmetryDesc.vtemp.x * pa1.x - JS.SymmetryDesc.vtemp.y * pa1.y - JS.SymmetryDesc.vtemp.z * pa1.z;
plane = JU.P4.new4(JS.SymmetryDesc.vtemp.x, JS.SymmetryDesc.vtemp.y, JS.SymmetryDesc.vtemp.z, w);
margin = (Math.abs(w) < 0.01 && JS.SymmetryDesc.vtemp.x * JS.SymmetryDesc.vtemp.y > 0.4 ? 1.30 : 1.05);
isRotation = false;
haveInversion = false;
isMirrorPlane = true;
break;
default:
haveInversion = false;
break;
}
if (f != 0) {
JS.SymmetryDesc.vtemp.sub2(pta00, pa1);
JS.SymmetryDesc.vtemp.add(pt0);
JS.SymmetryDesc.vtemp.sub(pa1);
JS.SymmetryDesc.vtemp.sub(d);
JS.SymmetryDesc.vtemp.scale(f);
pa1.add(JS.SymmetryDesc.vtemp);
ipt =  new JU.P3();
ipt.scaleAdd2(0.5, d, pa1);
ptinv =  new JU.P3();
ptinv.scaleAdd2(-2, ipt, pt0);
ptinv.scale(-1);
}} else if (trans != null) {
JS.SymmetryDesc.ptemp.setT(trans);
uc.toFractional(JS.SymmetryDesc.ptemp, false);
ftrans.setT(JS.SymmetryDesc.ptemp);
uc.toCartesian(JS.SymmetryDesc.ptemp, false);
trans.setT(JS.SymmetryDesc.ptemp);
}var ang = ang1;
JS.SymmetryDesc.approx0(ax1);
if (isRotation) {
var ptr =  new JU.P3();
JS.SymmetryDesc.vtemp.setT(ax1);
var ang2 = ang1;
var p0;
if (haveInversion) {
ptr.setT(ptinv);
p0 = ptinv;
} else if (pitch1 == 0) {
p0 = pt0;
ptr.setT(pa1);
} else {
p0 = pt0;
ptr.scaleAdd2(0.5, JS.SymmetryDesc.vtemp, pa1);
}JS.SymmetryDesc.ptemp.add2(pa1, JS.SymmetryDesc.vtemp);
ang2 = Math.round(JU.Measure.computeTorsion(pta00, pa1, JS.SymmetryDesc.ptemp, p0, true));
if (JS.SymmetryOperation.approx(ang2) != 0) {
ang1 = ang2;
if (ang1 < 0) ang1 = 360 + ang1;
}}var info1 = null;
var type = null;
var glideType = String.fromCharCode(0);
var order = op.getOpOrder();
var isccw = op.getOpIsCCW();
if (bsInfo.get(2) || bsInfo.get(15)) {
info1 = type = "identity";
if (isInversionOnly) {
JS.SymmetryDesc.ptemp.setT(ipt);
uc.toFractional(JS.SymmetryDesc.ptemp, false);
info1 = "Ci: " + JS.SymmetryDesc.strCoord(op, JS.SymmetryDesc.ptemp, op.isBio);
type = "inversion center";
} else if (isRotation) {
var screwtype = (isccw == null || haveInversion || pitch1 == 0 ? "" : isccw === Boolean.TRUE ? "(+)" : "(-)");
if (haveInversion) {
info1 = (Clazz.doubleToInt(360 / ang)) + "-bar" + screwtype + " axis";
} else if (pitch1 != 0) {
JS.SymmetryDesc.ptemp.setT(ax1);
uc.toFractional(JS.SymmetryDesc.ptemp, false);
info1 = (Clazz.doubleToInt(360 / ang)) + screwtype + " (" + JS.SymmetryDesc.strCoord(op, JS.SymmetryDesc.ptemp, op.isBio) + ") screw axis";
} else {
info1 = "C" + (Clazz.doubleToInt(360 / ang)) + screwtype + " axis";
}type = info1;
} else if (trans != null) {
var s = " " + JS.SymmetryDesc.strCoord(op, ftrans, op.isBio);
if (isTranslation) {
type = info1 = "translation";
info1 += ":" + s;
} else if (isMirrorPlane) {
if (isSpaceGroup) {
JS.SymmetryDesc.fixGlideTrans(ftrans);
trans.setT(ftrans);
uc.toCartesian(trans, true);
}s = " " + JS.SymmetryDesc.strCoord(op, ftrans, op.isBio);
glideType = JS.SymmetryOperation.getGlideFromTrans(ftrans, ax1);
type = info1 = glideType + "-glide plane";
info1 += "|translation:" + s;
}} else if (isMirrorPlane) {
type = info1 = "mirror plane";
}if (haveInversion && !isInversionOnly) {
JS.SymmetryDesc.ptemp.setT(ipt);
uc.toFractional(JS.SymmetryDesc.ptemp, false);
info1 += "|at " + JS.SymmetryDesc.strCoord(op, JS.SymmetryDesc.ptemp, op.isBio);
}if (isTimeReversed) {
info1 += "|time-reversed";
type += " (time-reversed)";
}}var isOK = true;
var isScrew = (isRotation && !haveInversion && pitch1 != 0);
if (!isScrew) {
isOK = JS.SymmetryDesc.checkHandedness(uc, ax1);
if (!isOK) {
ang1 = -ang1;
if (ang1 < 0) ang1 = 360 + ang1;
ax1.scale(-1);
}}var ignore = false;
var cmds = null;
if (id != null && bsInfo.get(3)) {
if (op.getOpType() == 0 || isSpaceGroupAll && op.isIrrelevant) {
System.out.println("!!SD irrelevent " + op.getOpTitle() + op.getOpPoint());
cmds = "";
} else {
var opType = null;
var drawid = "\ndraw ID " + id + "_";
var draw1 =  new JU.SB();
draw1.append(drawid).append("* delete");
if (!isSpaceGroup) {
JS.SymmetryDesc.drawLine(draw1, drawid + "frame1X", 0.15, pta00, JS.SymmetryDesc.pta01, "red");
JS.SymmetryDesc.drawLine(draw1, drawid + "frame1Y", 0.15, pta00, JS.SymmetryDesc.pta02, "green");
JS.SymmetryDesc.drawLine(draw1, drawid + "frame1Z", 0.15, pta00, pta03, "blue");
}var color;
var planeCenter = null;
var nPC = 0;
var isSpecial = (pta00.distance(pt0) < 0.2);
var title = (isSpaceGroup ? "<hover>" + id + ": " + op.xyz + "\n" + info1 + "</hover>" : null);
if (isRotation) {
color = "red";
ang = ang1;
var scale = 1;
JS.SymmetryDesc.vtemp.setT(ax1);
if (haveInversion) {
color = "purple";
}var wp = "";
if (isSpaceGroup) {
pa1.setT(op.getOpPoint());
uc.toCartesian(pa1, false);
}var ptr =  new JU.P3();
if (pitch1 != 0 && !haveInversion) {
opType = drawid + "screw";
color = (isccw === Boolean.TRUE ? "orange" : isccw === Boolean.FALSE ? "yellow" : order == 4 ? "lightgray" : "grey");
if (!isSpaceGroup) {
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine1", 0.1, pta00, pa1, "red");
JS.SymmetryDesc.ptemp.add2(pa1, JS.SymmetryDesc.vtemp);
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine2", 0.1, pt0, JS.SymmetryDesc.ptemp, "red");
ptr.scaleAdd2(0.5, JS.SymmetryDesc.vtemp, pa1);
}} else {
ptr.setT(pa1);
if (!isOK) {
if (!isSpecial && !isSpaceGroup) pa1.sub2(pa1, JS.SymmetryDesc.vtemp);
}if (haveInversion) {
opType = drawid + "rotinv";
if (isSpaceGroup) {
JS.SymmetryDesc.vtemp.normalize();
if (isccw === Boolean.TRUE) {
JS.SymmetryDesc.vtemp.scale(-1);
}} else {
if (pitch1 == 0) {
ptr.setT(ipt);
JS.SymmetryDesc.vtemp.scale(3 * scaleFactor);
if (isSpecial) {
JS.SymmetryDesc.ptemp.scaleAdd2(0.25, JS.SymmetryDesc.vtemp, pa1);
pa1.scaleAdd2(-0.24, JS.SymmetryDesc.vtemp, pa1);
ptr.scaleAdd2(0.31, JS.SymmetryDesc.vtemp, ptr);
color = "cyan";
} else {
JS.SymmetryDesc.ptemp.scaleAdd2(-1, JS.SymmetryDesc.vtemp, pa1);
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine1", 0.1, pta00, ipt, "red");
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine2", 0.1, ptinv, ipt, "red");
}} else {
if (!isSpecial) {
scale = pta00.distance(ptr);
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine1", 0.1, pta00, ptr, "red");
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine2", 0.1, ptinv, ptr, "red");
}}}} else {
opType = drawid + "rot";
JS.SymmetryDesc.vtemp.scale(3 * scaleFactor);
if (isSpecial) {
} else {
if (!isSpaceGroup) {
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine1", 0.1, pta00, ptr, "red");
JS.SymmetryDesc.drawLine(draw1, drawid + "rotLine2", 0.1, pt0, ptr, "red");
}}ptr.setT(pa1);
if (pitch1 == 0 && isSpecial) ptr.scaleAdd2(0.25, JS.SymmetryDesc.vtemp, ptr);
}}if (!isSpaceGroup) {
if (ang > 180) {
ang = 180 - ang;
}JS.SymmetryDesc.ptemp.add2(ptr, JS.SymmetryDesc.vtemp);
draw1.append(drawid).append("rotRotArrow arrow width 0.1 scale " + JU.PT.escF(scale) + " arc ").append(JU.Escape.eP(ptr)).append(JU.Escape.eP(JS.SymmetryDesc.ptemp));
JS.SymmetryDesc.ptemp.setT(pta00);
if (JS.SymmetryDesc.ptemp.distance(pt0) < 0.1) JS.SymmetryDesc.ptemp.set(Math.random(), Math.random(), Math.random());
draw1.append(JU.Escape.eP(JS.SymmetryDesc.ptemp));
JS.SymmetryDesc.ptemp.set(0, ang - 5 * Math.signum(ang), 0);
draw1.append(JU.Escape.eP(JS.SymmetryDesc.ptemp)).append(" color red");
}var d;
ignore = ignore || !op.opIsLong && (isSpaceGroupAll && pitch1 > 0 && !haveInversion && op.getOpTrans().length() > (order == 2 ? 0.71 : order == 3 ? 0.578 : order == 4 ? 0.51 : 0.3));
if (ignore) {
System.out.println("SD ignoring " + op.getOpTrans().length() + " " + op.getOpTitle() + op.xyz);
}var p2 = null;
var wscale = 1;
if (pitch1 == 0 && !haveInversion) {
JS.SymmetryDesc.ptemp.scaleAdd2(0.5, JS.SymmetryDesc.vtemp, pa1);
pa1.scaleAdd2(isSpaceGroup ? -0.5 : -0.45, JS.SymmetryDesc.vtemp, pa1);
if (isSpaceGroupAll && (p2 = op.getOpPoint2()) != null) {
ptr.setT(p2);
uc.toCartesian(ptr, false);
ptr.scaleAdd2(-0.5, JS.SymmetryDesc.vtemp, ptr);
}if (isSpaceGroup) {
this.scaleByOrder(JS.SymmetryDesc.vtemp, order, isccw);
}} else if (isSpaceGroupAll && pitch1 != 0 && !haveInversion && (d = op.getOpTrans().length()) > 0.4) {
if (isccw === Boolean.TRUE) {
} else if (isccw == null) {
} else if (d == 0.5) {
ignore = true;
}} else if (isSpaceGroup && haveInversion) {
p2 = pt1;
JS.SymmetryDesc.vtemp.scale(-1.0);
this.scaleByOrder(JS.SymmetryDesc.vtemp, order, isccw);
wp = "80";
}if (pitch1 > 0 && !haveInversion) {
wp = "" + (90 - Clazz.floatToInt(JS.SymmetryDesc.vtemp.length() * wscale / pitch1 * 90));
}if (!ignore) {
var name = opType + order + "rotVector1";
JS.SymmetryDesc.drawVector(draw1, drawid, name, "vector", "0.1" + wp, pa1, JS.SymmetryDesc.vtemp, isTimeReversed ? "gray" : color, title);
if (p2 != null) {
JS.SymmetryDesc.drawVector(draw1, drawid, name + "b", "vector", "0.1" + wp, ptr, JS.SymmetryDesc.vtemp, isTimeReversed ? "gray" : color, title);
}}} else if (isMirrorPlane) {
JS.SymmetryDesc.ptemp.sub2(ptref, pta00);
if (!isSpaceGroup && pta00.distance(ptref) > 0.2) JS.SymmetryDesc.drawVector(draw1, drawid, "planeVector", "vector", "0.05", pta00, JS.SymmetryDesc.ptemp, isTimeReversed ? "gray" : "cyan", null);
opType = drawid + "plane";
var p = JU.P4.newPt(plane);
if (trans == null) {
color = "green";
} else {
opType = drawid + "glide";
switch ((glideType).charCodeAt(0)) {
case 103:
color = "gold";
break;
case 97:
color = "magenta";
break;
case 98:
color = "violet";
break;
case 99:
color = "blue";
break;
case 110:
color = "orange";
break;
case 100:
color = "grey";
break;
default:
color = "grey";
break;
}
p.w += 0.01;
if (!isSpaceGroup) {
JS.SymmetryDesc.drawFrameLine("X", ptref, vt1, 0.15, JS.SymmetryDesc.ptemp, draw1, opType, "red");
JS.SymmetryDesc.drawFrameLine("Y", ptref, vt2, 0.15, JS.SymmetryDesc.ptemp, draw1, opType, "green");
JS.SymmetryDesc.drawFrameLine("Z", ptref, vt3, 0.15, JS.SymmetryDesc.ptemp, draw1, opType, "blue");
}}var v = this.modelSet.vwr.getTriangulator().intersectPlane(p, uc.getCanonicalCopy(margin, true), 3);
if (v != null) {
var isCoincident = (isSpaceGroup && op.isCoincident);
planeCenter =  new JU.P3();
var pt = draw1.length();
for (var i = v.size(); --i >= 0; ) {
var pts = v.get(i);
draw1.append(drawid).append(trans == null ? "m" : "g").append("planep").appendI(i).append(" ").append(JU.Escape.eP(pts[0])).append(JU.Escape.eP(pts[1]));
if (pts.length == 3) {
if (!isCoincident || (i % 2 == 0) != (trans == null)) {
draw1.append(JU.Escape.eP(pts[2]));
}} else {
planeCenter.add(pts[0]);
planeCenter.add(pts[1]);
nPC += 2;
}draw1.append(" color translucent ").append(color);
if (title != null) draw1.append(" ").append(JU.PT.esc(title));
}
}if (v == null || v.size() == 0) {
if (isSpaceGroupAll) {
ignore = true;
} else {
JS.SymmetryDesc.ptemp.add2(pa1, ax1);
draw1.append(drawid).append("planeCircle scale 2.0 circle ").append(JU.Escape.eP(pa1)).append(JU.Escape.eP(JS.SymmetryDesc.ptemp)).append(" color translucent ").append(color).append(" mesh fill");
if (title != null) draw1.append(" ").append(JU.PT.esc(title));
}}}if (haveInversion) {
opType = drawid + "inv";
if (isInversionOnly) {
draw1.append(drawid).append("invPoint diameter 0.4 ").append(JU.Escape.eP(ipt));
if (title != null) draw1.append(" ").append(JU.PT.esc(title));
JS.SymmetryDesc.ptemp.sub2(ptinv, pta00);
if (!isSpaceGroup) {
JS.SymmetryDesc.drawVector(draw1, drawid, "Arrow", "vector", "0.05", pta00, JS.SymmetryDesc.ptemp, isTimeReversed ? "gray" : "cyan", null);
}} else {
if (order == 4) {
draw1.append(drawid).append("RotPoint diameter 0.3 color red").append(JU.Escape.eP(ipt));
if (title != null) draw1.append(" ").append(JU.PT.esc(title));
}if (!isSpaceGroup) {
draw1.append(" color cyan");
if (!isSpecial) {
JS.SymmetryDesc.ptemp.sub2(pt0, ptinv);
JS.SymmetryDesc.drawVector(draw1, drawid, "Arrow", "vector", "0.05", ptinv, JS.SymmetryDesc.ptemp, isTimeReversed ? "gray" : "cyan", null);
}if (options != 1073742066) {
JS.SymmetryDesc.vtemp.setT(vt1);
JS.SymmetryDesc.vtemp.scale(-1);
JS.SymmetryDesc.drawFrameLine("X", ptinv, JS.SymmetryDesc.vtemp, 0.15, JS.SymmetryDesc.ptemp, draw1, opType, "red");
JS.SymmetryDesc.vtemp.setT(vt2);
JS.SymmetryDesc.vtemp.scale(-1);
JS.SymmetryDesc.drawFrameLine("Y", ptinv, JS.SymmetryDesc.vtemp, 0.15, JS.SymmetryDesc.ptemp, draw1, opType, "green");
JS.SymmetryDesc.vtemp.setT(vt3);
JS.SymmetryDesc.vtemp.scale(-1);
JS.SymmetryDesc.drawFrameLine("Z", ptinv, JS.SymmetryDesc.vtemp, 0.15, JS.SymmetryDesc.ptemp, draw1, opType, "blue");
}}}}if (trans != null) {
if (isMirrorPlane && isSpaceGroup) {
if (planeCenter != null) {
ptref = planeCenter;
ptref.scale(1 / nPC);
ptref.scaleAdd2(-0.5, trans, ptref);
}} else if (ptref == null) {
ptref = (isSpaceGroup ?  new JU.P3() : JU.P3.newP(pta00));
}if (ptref != null && !ignore) {
JS.SymmetryDesc.drawVector(draw1, drawid, "transVector", "vector", (isTranslationOnly ? "0.1" : "0.05"), ptref, trans, isTimeReversed && !haveInversion && !isMirrorPlane && !isRotation ? "darkGray" : "gold", title);
}}if (!isSpaceGroup) {
JS.SymmetryDesc.ptemp2.setT(pt0);
JS.SymmetryDesc.ptemp.sub2(pt1, pt0);
JS.SymmetryDesc.ptemp.scaleAdd2(0.9, JS.SymmetryDesc.ptemp, JS.SymmetryDesc.ptemp2);
JS.SymmetryDesc.drawLine(draw1, drawid + "frame2X", 0.2, JS.SymmetryDesc.ptemp2, JS.SymmetryDesc.ptemp, "red");
JS.SymmetryDesc.ptemp.sub2(pt2, pt0);
JS.SymmetryDesc.ptemp.scaleAdd2(0.9, JS.SymmetryDesc.ptemp, JS.SymmetryDesc.ptemp2);
JS.SymmetryDesc.drawLine(draw1, drawid + "frame2Y", 0.2, JS.SymmetryDesc.ptemp2, JS.SymmetryDesc.ptemp, "green");
JS.SymmetryDesc.ptemp.sub2(pt3, pt0);
JS.SymmetryDesc.ptemp.scaleAdd2(0.9, JS.SymmetryDesc.ptemp, JS.SymmetryDesc.ptemp2);
JS.SymmetryDesc.drawLine(draw1, drawid + "frame2Z", 0.2, JS.SymmetryDesc.ptemp2, JS.SymmetryDesc.ptemp, "purple");
draw1.append("\nsym_point = " + JU.Escape.eP(pta00));
draw1.append("\nvar p0 = " + JU.Escape.eP(JS.SymmetryDesc.ptemp2));
if (Clazz.instanceOf(pta00,"JM.Atom")) {
draw1.append("\nvar set2 = within(0.2,p0);if(!set2){set2 = within(0.2,p0.uxyz.xyz)}");
draw1.append("\n set2 &= {_" + (pta00).getElementSymbol() + "}");
} else {
draw1.append("\nvar set2 = p0.uxyz");
}draw1.append("\nsym_target = set2;if (set2) {");
if (!isSpecial && options != 1073742066 && ptTarget == null && !haveTranslation) {
draw1.append(drawid).append("offsetFrameX diameter 0.20 @{set2.xyz} @{set2.xyz + ").append(JU.Escape.eP(vt1)).append("*0.9} color red");
draw1.append(drawid).append("offsetFrameY diameter 0.20 @{set2.xyz} @{set2.xyz + ").append(JU.Escape.eP(vt2)).append("*0.9} color green");
draw1.append(drawid).append("offsetFrameZ diameter 0.20 @{set2.xyz} @{set2.xyz + ").append(JU.Escape.eP(vt3)).append("*0.9} color purple");
}draw1.append("\n}\n");
}cmds = draw1.toString();
if (JU.Logger.debugging) JU.Logger.info(cmds);
draw1 = null;
drawid = null;
}}if (trans == null) ftrans = null;
if (isScrew) {
trans = JU.V3.newV(ax1);
JS.SymmetryDesc.ptemp.setT(trans);
uc.toFractional(JS.SymmetryDesc.ptemp, false);
ftrans = JU.V3.newV(JS.SymmetryDesc.ptemp);
}if (isMirrorPlane) {
ang1 = 0;
}if (haveInversion) {
if (isInversionOnly) {
pa1 = null;
ax1 = null;
trans = null;
ftrans = null;
}} else if (isTranslation) {
pa1 = null;
ax1 = null;
}if (ax1 != null) ax1.normalize();
var xyzNew = null;
if (bsInfo.get(0) || bsInfo.get(17)) {
xyzNew = (op.isBio ? m2.toString() : op.modDim > 0 ? op.xyzOriginal : JS.SymmetryOperation.getXYZFromMatrix(m2, false, false, false));
if (isMagnetic) xyzNew = op.fixMagneticXYZ(m2, xyzNew, true);
}var ret =  new Array(20);
for (var i = bsInfo.nextSetBit(0); i >= 0; i = bsInfo.nextSetBit(i + 1)) {
switch (i) {
case 0:
ret[i] = xyzNew;
break;
case 19:
if (ptFrom != null && ptTarget == null && !op.isBio && op.modDim == 0) {
var xyzN;
JS.SymmetryDesc.pta02.setT(ptFrom);
uc.toFractional(JS.SymmetryDesc.pta02, true);
m2.rotTrans(JS.SymmetryDesc.pta02);
JS.SymmetryDesc.ptemp.setT(JS.SymmetryDesc.pta02);
uc.unitize(JS.SymmetryDesc.pta02);
JS.SymmetryDesc.vtrans.sub2(JS.SymmetryDesc.pta02, JS.SymmetryDesc.ptemp);
m2 = JU.M4.newM4(op);
m2.add(JS.SymmetryDesc.vtrans);
xyzN = JS.SymmetryOperation.getXYZFromMatrix(m2, false, false, false);
if (isMagnetic) xyzN = op.fixMagneticXYZ(m2, xyzN, true);
ret[i] = xyzN;
}break;
case 1:
ret[i] = op.xyzOriginal;
break;
case 2:
ret[i] = info1;
break;
case 3:
ret[i] = cmds;
break;
case 4:
ret[i] = JS.SymmetryDesc.approx0(ftrans);
break;
case 5:
ret[i] = JS.SymmetryDesc.approx0(trans);
break;
case 6:
ret[i] = JS.SymmetryDesc.approx0(ipt);
break;
case 7:
ret[i] = JS.SymmetryDesc.approx0(pa1 != null && bsInfo.get(22) ? pta00 : pa1);
break;
case 8:
ret[i] = (plane == null ? JS.SymmetryDesc.approx0(ax1) : null);
break;
case 9:
ret[i] = (ang1 != 0 ? Integer.$valueOf(ang1) : null);
break;
case 10:
ret[i] = m2;
break;
case 11:
ret[i] = (JS.SymmetryDesc.vtrans.lengthSquared() > 0 ? JS.SymmetryDesc.vtrans : null);
break;
case 12:
ret[i] = op.getCentering();
break;
case 13:
ret[i] = Integer.$valueOf(op.timeReversal);
break;
case 14:
if (plane != null && bsInfo.get(22)) {
var d = JU.Measure.distanceToPlane(plane, pta00);
plane.w -= d;
}ret[i] = plane;
break;
case 15:
ret[i] = type;
break;
case 16:
ret[i] = Integer.$valueOf(op.number);
break;
case 17:
var cift = null;
if (!op.isBio && !xyzNew.equals(op.xyzOriginal)) {
if (op.number > 0) {
var orig = JS.SymmetryOperation.getMatrixFromXYZ(op.xyzOriginal, null, false);
orig.sub(m2);
cift =  new JU.P3();
orig.getTranslation(cift);
}}var cifi = (op.number < 0 ? 0 : op.number);
ret[i] = cifi + (cift == null ? " [0 0 0]" : " [" + Clazz.floatToInt(-cift.x) + " " + Clazz.floatToInt(-cift.y) + " " + Clazz.floatToInt(-cift.z) + "]");
break;
case 18:
ret[i] = op.xyzCanonical;
break;
}
}
return ret;
}, "JS.SymmetryOperation,J.api.SymmetryInterface,JU.P3,JU.P3,~S,~N,~N,~B,JU.BS,~B,~B");
c$.fixGlideTrans = Clazz.defineMethod(c$, "fixGlideTrans", 
function(ftrans){
ftrans.x = JS.SymmetryDesc.fixGlideX(ftrans.x);
ftrans.y = JS.SymmetryDesc.fixGlideX(ftrans.y);
ftrans.z = JS.SymmetryDesc.fixGlideX(ftrans.z);
}, "JU.V3");
c$.fixGlideX = Clazz.defineMethod(c$, "fixGlideX", 
function(x){
switch (Math.round(x * 12.001)) {
case 9:
return -0.25;
case -9:
return 0.25;
default:
return x;
}
}, "~N");
c$.checkHandedness = Clazz.defineMethod(c$, "checkHandedness", 
function(uc, ax1){
var a;
var b;
var c;
JS.SymmetryDesc.ptemp.set(1, 0, 0);
uc.toCartesian(JS.SymmetryDesc.ptemp, false);
a = JS.SymmetryDesc.approx0d(JS.SymmetryDesc.ptemp.dot(ax1));
JS.SymmetryDesc.ptemp.set(0, 1, 0);
uc.toCartesian(JS.SymmetryDesc.ptemp, false);
b = JS.SymmetryDesc.approx0d(JS.SymmetryDesc.ptemp.dot(ax1));
JS.SymmetryDesc.ptemp.set(0, 0, 1);
uc.toCartesian(JS.SymmetryDesc.ptemp, false);
c = JS.SymmetryDesc.approx0d(JS.SymmetryDesc.ptemp.dot(ax1));
return (a == 0 ? (b == 0 ? c > 0 : b > 0) : c == 0 ? a > 0 : (b == 0 ? c > 0 : a * b * c > 0));
}, "J.api.SymmetryInterface,JU.P3");
Clazz.defineMethod(c$, "scaleByOrder", 
function(v, order, isccw){
v.scale(1 + (0.3 / order) + (isccw == null ? 0 : isccw === Boolean.TRUE ? 0.02 : -0.02));
}, "JU.V3,~N,Boolean");
c$.drawLine = Clazz.defineMethod(c$, "drawLine", 
function(s, id, diameter, pt0, pt1, color){
s.append(id).append(" diameter ").appendD(diameter).append(JU.Escape.eP(pt0)).append(JU.Escape.eP(pt1)).append(" color ").append(color);
}, "JU.SB,~S,~N,JU.P3,JU.P3,~S");
c$.drawFrameLine = Clazz.defineMethod(c$, "drawFrameLine", 
function(xyz, pt, v, width, ptemp, draw1, key, color){
ptemp.setT(pt);
ptemp.add(v);
JS.SymmetryDesc.drawLine(draw1, key + "Pt" + xyz, width, pt, ptemp, "translucent " + color);
}, "~S,JU.P3,JU.V3,~N,JU.P3,JU.SB,~S,~S");
c$.drawVector = Clazz.defineMethod(c$, "drawVector", 
function(draw1, drawid, label, type, d, pt1, v, color, title){
if (type.equals("vline")) {
JS.SymmetryDesc.ptemp2.add2(pt1, v);
type = "";
v = JS.SymmetryDesc.ptemp2;
}d += " ";
draw1.append(drawid).append(label).append(" diameter ").append(d).append(type).append(JU.Escape.eP(pt1)).append(JU.Escape.eP(v)).append(" color ").append(color);
if (title != null) draw1.append(" \"" + title + "\"");
}, "JU.SB,~S,~S,~S,~S,JU.T3,JU.T3,~S,~S");
c$.setFractional = Clazz.defineMethod(c$, "setFractional", 
function(uc, pt00, pt01, offset){
pt01.setT(pt00);
if (offset != null) uc.toUnitCell(pt01, offset);
uc.toFractional(pt01, false);
}, "J.api.SymmetryInterface,JU.T3,JU.P3,JU.P3");
c$.rotTransCart = Clazz.defineMethod(c$, "rotTransCart", 
function(op, uc, pt00, vtrans){
var p0 = JU.P3.newP(pt00);
uc.toFractional(p0, false);
op.rotTrans(p0);
p0.add(vtrans);
uc.toCartesian(p0, false);
return p0;
}, "JS.SymmetryOperation,J.api.SymmetryInterface,JU.P3,JU.V3");
c$.strCoord = Clazz.defineMethod(c$, "strCoord", 
function(op, p, isBio){
JS.SymmetryDesc.approx0(p);
return (isBio ? p.x + " " + p.y + " " + p.z : op.fcoord2(p));
}, "JS.SymmetryOperation,JU.T3,~B");
c$.approx0 = Clazz.defineMethod(c$, "approx0", 
function(pt){
if (pt != null) {
pt.x = JS.SymmetryDesc.approx0d(pt.x);
pt.y = JS.SymmetryDesc.approx0d(pt.y);
pt.z = JS.SymmetryDesc.approx0d(pt.z);
}return pt;
}, "JU.T3");
c$.approx0d = Clazz.defineMethod(c$, "approx0d", 
function(x){
return (Math.abs(x) < 0.0001 ? 0 : x);
}, "~N");
Clazz.defineMethod(c$, "getSymmetryInfo", 
function(iModel, iatom, uc, xyz, op, translation, pt, pt2, id, type, scaleFactor, nth, options, isSpaceGroup){
var returnType = 0;
var nullRet = JS.SymmetryDesc.nullReturn(type);
switch (type) {
case 1073741994:
return "" + uc.getLatticeType();
case 1073742001:
returnType = 1825200146;
break;
case 135176:
returnType = 135176;
break;
case 1275068418:
returnType = JS.SymmetryDesc.getType(id);
switch (returnType) {
case 1153433601:
case 1073741961:
case 1073742001:
case 134217751:
case 1086326789:
case 36868:
type = returnType;
break;
default:
returnType = JS.SymmetryDesc.getKeyType(id);
break;
}
break;
}
var bsInfo = JS.SymmetryDesc.getInfoBS(returnType);
var isSpaceGroupAll = (nth == -2);
var iop = op;
var offset = (options == 1073742066 && (type == 1153433601 || type == 134217751) ? pt2 : null);
if (offset != null) pt2 = null;
var info = null;
var xyzOriginal = null;
var ops = null;
if (pt2 == null) {
if (xyz == null) {
ops = (isSpaceGroupAll ? uc.getAdditionalOperations() : uc.getSymmetryOperations());
if (ops == null || Math.abs(op) > ops.length) return nullRet;
if (op == 0) return nullRet;
iop = Math.abs(op) - 1;
xyz = (translation == null ? ops[iop].xyz : ops[iop].getxyzTrans(translation));
xyzOriginal = ops[iop].xyzOriginal;
} else {
iop = op = 0;
}var symTemp =  new JS.Symmetry();
symTemp.setSpaceGroup(false);
var isBio = (uc != null && uc.isBio());
var i = (isBio ? symTemp.addBioMoleculeOperation((uc.getSpaceGroup()).finalOperations[iop], op < 0) : symTemp.addSpaceGroupOperation((op < 0 ? "!" : "=") + xyz, Math.abs(op)));
if (i < 0) return nullRet;
var opTemp = symTemp.getSpaceGroupOperation(i);
if (isSpaceGroup) opTemp.isCoincident = ops[iop].isCoincident;
if (xyzOriginal != null) opTemp.xyzOriginal = xyzOriginal;
opTemp.number = op;
if (!isBio) opTemp.getCentering();
if (pt == null && iatom >= 0) pt = this.modelSet.at[iatom];
if (type == 134217751 || type == 1153433601) {
if (isBio || pt == null) return nullRet;
symTemp.setUnitCell(uc);
JS.SymmetryDesc.ptemp.setT(pt);
uc.toFractional(JS.SymmetryDesc.ptemp, false);
if (Float.isNaN(JS.SymmetryDesc.ptemp.x)) return nullRet;
var sympt =  new JU.P3();
symTemp.newSpaceGroupPoint(JS.SymmetryDesc.ptemp, i, null, 0, 0, 0, sympt);
if (options == 1073742066) {
uc.unitize(sympt);
sympt.add(offset);
}symTemp.toCartesian(sympt, false);
var ret = sympt;
return (type == 1153433601 ? this.getAtom(uc, iModel, iatom, ret) : ret);
}info = this.createInfoArray(opTemp, uc, pt, null, (id == null ? "sym" : id), scaleFactor, options, (translation != null), bsInfo, isSpaceGroup, isSpaceGroupAll);
if (type == 1275068418 && id != null) {
returnType = JS.SymmetryDesc.getKeyType(id);
}} else {
var stype = "info";
var asString = false;
switch (type) {
case 1275068418:
returnType = JS.SymmetryDesc.getKeyType(id);
id = stype = null;
if (nth == 0) nth = -1;
break;
case 1073742001:
id = stype = null;
if (nth == 0) nth = -1;
asString = true;
bsInfo.set(21);
bsInfo.set(0);
bsInfo.set(19);
break;
case 135176:
if (id == null) id = (isSpaceGroup ? "sg" : "sym");
stype = "all";
asString = true;
break;
case 1153433601:
id = stype = null;
default:
if (nth == 0) nth = 1;
}
var ret1 = this.getSymopInfoForPoints(uc, iModel, op, translation, pt, pt2, id, stype, scaleFactor, nth, options, bsInfo);
if (asString) {
return ret1;
}if ((typeof(ret1)=='string')) return nullRet;
info = ret1;
if (type == 1153433601) {
if (!(Clazz.instanceOf(pt,"JM.Atom")) && !(Clazz.instanceOf(pt2,"JM.Atom"))) iatom = -1;
return (info == null ? nullRet : this.getAtom(uc, iModel, iatom, info[7]));
}}if (info == null) return nullRet;
var isList = (info.length > 0 && Clazz.instanceOf(info[0],Array));
if (nth < 0 && op <= 0 && xyz == null && (type == 1275068418 || isList)) {
if (type == 1275068418 && info.length > 0 && !(Clazz.instanceOf(info[0],Array))) info =  Clazz.newArray(-1, [info]);
var lst =  new JU.Lst();
for (var i = 0; i < info.length; i++) lst.addLast(JS.SymmetryDesc.getInfo(info[i], returnType < 0 ? returnType : type));

return lst;
} else if (returnType < 0 && (nth >= 0 || op > 0 || xyz != null)) {
type = returnType;
}if (nth > 0 && isList) info = info[0];
if (type == 135176 && isSpaceGroup && nth == -2) type = 134217764;
return JS.SymmetryDesc.getInfo(info, type);
}, "~N,~N,J.api.SymmetryInterface,~S,~N,JU.P3,JU.P3,JU.P3,~S,~N,~N,~N,~N,~B");
Clazz.defineMethod(c$, "getAtom", 
function(uc, iModel, iAtom, sympt){
var bsElement = null;
if (iAtom >= 0) this.modelSet.getAtomBitsMDa(1094715402, Integer.$valueOf(this.modelSet.at[iAtom].getElementNumber()), bsElement =  new JU.BS());
var bsResult =  new JU.BS();
this.modelSet.getAtomsWithin(0.02, sympt, bsResult, iModel);
if (bsElement != null) bsResult.and(bsElement);
if (bsResult.isEmpty()) {
sympt = JU.P3.newP(sympt);
uc.toUnitCell(sympt, null);
uc.toCartesian(sympt, false);
this.modelSet.getAtomsWithin(0.02, sympt, bsResult, iModel);
if (bsElement != null) bsResult.and(bsElement);
}return bsResult;
}, "J.api.SymmetryInterface,~N,~N,JU.T3");
Clazz.defineMethod(c$, "getSymopInfoForPoints", 
function(sym, modelIndex, symOp, translation, pt1, pt2, drawID, stype, scaleFactor, nth, options, bsInfo){
var asString = (bsInfo.get(21) || bsInfo.get(3) && bsInfo.cardinality() == 3);
bsInfo.clear(21);
var ret = (asString ? "" : null);
var sginfo = this.getSpaceGroupInfo(sym, modelIndex, null, symOp, pt1, pt2, drawID, scaleFactor, nth, false, true, options, null, bsInfo);
if (sginfo == null) return ret;
var infolist = sginfo.get("operations");
if (infolist == null) return ret;
var sb = (asString ?  new JU.SB() : null);
symOp--;
var isAll = (!asString && symOp < 0);
var strOperations = sginfo.get("symmetryInfo");
var labelOnly = "label".equals(stype);
var n = 0;
for (var i = 0; i < infolist.length; i++) {
if (infolist[i] == null || symOp >= 0 && symOp != i) continue;
if (!asString) {
if (!isAll) return infolist[i];
infolist[n++] = infolist[i];
continue;
}if (drawID != null) return (infolist[i][3]) + "\nprint " + JU.PT.esc(strOperations);
if (sb.length() > 0) sb.appendC('\n');
if (!labelOnly) {
if (symOp < 0) sb.appendI(i + 1).appendC('\t');
sb.append(infolist[i][0]).appendC('\t');
}sb.append(infolist[i][2]);
}
if (!asString) {
var a =  new Array(n);
for (var i = 0; i < n; i++) a[i] = infolist[i];

return a;
}if (sb.length() == 0) return (drawID != null ? "draw " + drawID + "* delete" : ret);
return sb.toString();
}, "J.api.SymmetryInterface,~N,~N,JU.P3,JU.P3,JU.P3,~S,~S,~N,~N,~N,JU.BS");
Clazz.defineMethod(c$, "getSymopInfo", 
function(iAtom, xyz, op, translation, pt, pt2, id, type, scaleFactor, nth, options, opList){
if (type == 0) type = JS.SymmetryDesc.getType(id);
var ret = (type == 1153433601 ?  new JU.BS() : "");
var iModel = (iAtom >= 0 ? this.modelSet.at[iAtom].mi : this.modelSet.vwr.am.cmi);
if (iModel < 0) return ret;
var uc = this.modelSet.am[iModel].biosymmetry;
if (uc == null && (uc = this.modelSet.getUnitCell(iModel)) == null) {
uc =  new JS.Symmetry().setUnitCell(null, false, NaN);
}if (type != 135176 || op != 2147483647 && opList == null) {
return this.getSymmetryInfo(iModel, iAtom, uc, xyz, op, translation, pt, pt2, id, type, scaleFactor, nth, options, false);
}if (uc == null) return ret;
var isSpaceGroup = (xyz == null && nth < 0 && opList == null);
var s = "";
var ops = (isSpaceGroup && nth == -2 ? uc.getAdditionalOperations() : uc.getSymmetryOperations());
if (ops != null) {
if (id == null) id = "sg";
var n = ops.length;
if (pt != null && pt2 == null || opList != null) {
if (opList == null) opList = uc.getInvariantSymops(pt, null);
n = opList.length;
for (var i = 0; i < n; i++) {
if (nth > 0 && nth != i + 1) continue;
op = opList[i];
s += this.getSymmetryInfo(iModel, iAtom, uc, xyz, op, translation, pt, pt2, id + op, 135176, scaleFactor, nth, options, pt == null);
}
} else {
for (op = 1; op <= n; op++) s += this.getSymmetryInfo(iModel, iAtom, uc, xyz, op, translation, pt, pt2, id + op, 135176, scaleFactor, nth, options, true);

}}return s;
}, "~N,~S,~N,JU.P3,JU.P3,JU.P3,~S,~N,~N,~N,~N,~A");
Clazz.defineMethod(c$, "getSpaceGroupInfo", 
function(sym, modelIndex, sgName, symOp, pt1, pt2, drawID, scaleFactor, nth, isFull, isForModel, options, cellInfo, bsInfo){
if (bsInfo == null) {
bsInfo =  new JU.BS();
bsInfo.setBits(0, JS.SymmetryDesc.keys.length);
bsInfo.clear(19);
}var matrixOnly = (bsInfo.cardinality() == 1 && bsInfo.get(10));
var info = null;
var isStandard = (!matrixOnly && pt1 == null && drawID == null && nth <= 0 && bsInfo.cardinality() >= JS.SymmetryDesc.keys.length);
var isBio = false;
var sgNote = null;
var haveName = (sgName != null && sgName.length > 0);
var haveRawName = (haveName && sgName.indexOf("[--]") >= 0);
if (isForModel || !haveName) {
var saveModelInfo = (isStandard && symOp == 0);
if (matrixOnly) {
cellInfo = sym;
} else {
if (modelIndex < 0) modelIndex = (Clazz.instanceOf(pt1,"JM.Atom") ? (pt1).mi : this.modelSet.vwr.am.cmi);
if (modelIndex < 0) sgNote = "no single current model";
 else if (cellInfo == null && !(isBio = (cellInfo = this.modelSet.am[modelIndex].biosymmetry) != null) && (cellInfo = this.modelSet.getUnitCell(modelIndex)) == null) sgNote = "not applicable";
if (sgNote != null) {
info =  new java.util.Hashtable();
info.put("spaceGroupInfo", "");
info.put("spaceGroupNote", sgNote);
info.put("symmetryInfo", "");
} else if (isStandard) {
info = this.modelSet.getInfo(modelIndex, "spaceGroupInfo");
}if (info != null) return info;
sgName = cellInfo.getSpaceGroupName();
}info =  new java.util.Hashtable();
var ops = cellInfo.getSymmetryOperations();
var sg = (isBio ? (cellInfo).spaceGroup : null);
var slist = (haveRawName ? "" : null);
var opCount = 0;
if (ops != null) {
if (!matrixOnly) {
if (isBio) sym.setSpaceGroupTo(JS.SpaceGroup.getNull(false, false, false));
 else sym.setSpaceGroup(false);
}if (ops[0].timeReversal != 0) (sym.getSpaceGroupOperation(0)).timeReversal = 1;
var infolist =  new Array(ops.length);
var sops = "";
for (var i = 0, nop = 0; i < ops.length && nop != nth; i++) {
var op = ops[i];
var xyzOriginal = op.xyzOriginal;
var iop;
if (matrixOnly) {
iop = i;
} else {
var isNewIncomm = (i == 0 && op.xyz.indexOf("x4") >= 0);
iop = (!isNewIncomm && sym.getSpaceGroupOperation(i) != null ? i : isBio ? sym.addBioMoleculeOperation(sg.finalOperations[i], false) : sym.addSpaceGroupOperation("=" + op.xyz, i + 1));
if (iop < 0) continue;
op = sym.getSpaceGroupOperation(i);
if (op == null) continue;
op.xyzOriginal = xyzOriginal;
}if (op.timeReversal != 0 || op.modDim > 0) isStandard = false;
if (slist != null) slist += ";" + op.xyz;
var ret = (symOp > 0 && symOp - 1 != iop ? null : this.createInfoArray(op, cellInfo, pt1, pt2, drawID, scaleFactor, options, false, bsInfo, false, false));
if (ret != null) {
if (nth > 0 && ++nop != nth) continue;
infolist[i] = ret;
if (!matrixOnly) sops += "\n" + (i + 1) + "\t" + ret[bsInfo.get(19) ? 19 : 0] + "\t  " + ret[2];
opCount++;
if (symOp > 0) break;
}}
info.put("operations", infolist);
if (!matrixOnly) info.put("symmetryInfo", (sops.length == 0 ? "" : sops.substring(1)));
}if (matrixOnly) {
return info;
}sgNote = (opCount == 0 ? "\n no symmetry operations" : nth <= 0 && symOp <= 0 ? "\n" + opCount + " symmetry operation" + (opCount == 1 ? ":\n" : "s:\n") : "");
if (slist != null) sgName = slist.substring(slist.indexOf(";") + 1);
if (saveModelInfo) this.modelSet.setInfo(modelIndex, "spaceGroupInfo", info);
} else {
info =  new java.util.Hashtable();
}info.put("spaceGroupName", sgName);
info.put("spaceGroupNote", sgNote == null ? "" : sgNote);
var data;
if (isBio) {
data = sgName;
} else {
if (haveName && !haveRawName) sym.setSpaceGroupName(sgName);
data = sym.getSpaceGroupInfoObj(sgName, (cellInfo == null ? null : cellInfo.getUnitCellParams()), isFull, !isForModel);
if (data == null || data.equals("?")) {
data = "?";
info.put("spaceGroupNote", "could not identify space group from name: " + sgName + "\nformat: show spacegroup \"2\" or \"P 2c\" " + "or \"C m m m\" or \"x, y, z;-x ,-y, -z\"");
}}info.put("spaceGroupInfo", data);
return info;
}, "J.api.SymmetryInterface,~N,~S,~N,JU.P3,JU.P3,~S,~N,~N,~B,~B,~N,J.api.SymmetryInterface,JU.BS");
Clazz.defineMethod(c$, "getTransform", 
function(uc, ops, fracA, fracB, best){
if (JS.SymmetryDesc.pta01 == null) {
JS.SymmetryDesc.pta01 =  new JU.P3();
JS.SymmetryDesc.pta02 =  new JU.P3();
JS.SymmetryDesc.ptemp =  new JU.P3();
JS.SymmetryDesc.vtrans =  new JU.V3();
}JS.SymmetryDesc.pta02.setT(fracB);
JS.SymmetryDesc.vtrans.setT(JS.SymmetryDesc.pta02);
uc.unitize(JS.SymmetryDesc.pta02);
var dmin = 3.4028235E38;
var imin = -1;
for (var i = 0, n = ops.length; i < n; i++) {
var op = ops[i];
JS.SymmetryDesc.pta01.setT(fracA);
op.rotTrans(JS.SymmetryDesc.pta01);
JS.SymmetryDesc.ptemp.setT(JS.SymmetryDesc.pta01);
uc.unitize(JS.SymmetryDesc.pta01);
var d = JS.SymmetryDesc.pta01.distanceSquared(JS.SymmetryDesc.pta02);
if (d < 1.96E-6) {
JS.SymmetryDesc.vtrans.sub(JS.SymmetryDesc.ptemp);
JS.SymmetryOperation.normalize12ths(JS.SymmetryDesc.vtrans);
var m2 = JU.M4.newM4(op);
m2.add(JS.SymmetryDesc.vtrans);
JS.SymmetryDesc.pta01.setT(fracA);
m2.rotTrans(JS.SymmetryDesc.pta01);
uc.unitize(JS.SymmetryDesc.pta01);
d = JS.SymmetryDesc.pta01.distanceSquared(JS.SymmetryDesc.pta02);
if (d >= 1.96E-6) {
continue;
}return m2;
}if (d < dmin) {
dmin = d;
imin = i;
}}
if (best) {
var op = ops[imin];
JS.SymmetryDesc.pta01.setT(fracA);
op.rotTrans(JS.SymmetryDesc.pta01);
uc.unitize(JS.SymmetryDesc.pta01);
System.err.println("" + imin + " " + JS.SymmetryDesc.pta01.distance(JS.SymmetryDesc.pta02) + " " + JS.SymmetryDesc.pta01 + " " + JS.SymmetryDesc.pta02 + " " + JU.V3.newVsub(JS.SymmetryDesc.pta02, JS.SymmetryDesc.pta01));
}return null;
}, "JS.UnitCell,~A,JU.P3,JU.P3,~B");
c$.keys =  Clazz.newArray(-1, ["xyz", "xyzOriginal", "label", null, "fractionalTranslation", "cartesianTranslation", "inversionCenter", null, "axisVector", "rotationAngle", "matrix", "unitTranslation", "centeringVector", "timeReversal", "plane", "_type", "id", "cif2", "xyzCanonical", "xyzNormalized"]);
c$.vtemp =  new JU.V3();
c$.ptemp =  new JU.P3();
c$.ptemp2 =  new JU.P3();
c$.pta01 =  new JU.P3();
c$.pta02 =  new JU.P3();
c$.vtrans =  new JU.V3();
});
;//5.0.1-v2 Mon Feb 05 08:46:23 CST 2024
