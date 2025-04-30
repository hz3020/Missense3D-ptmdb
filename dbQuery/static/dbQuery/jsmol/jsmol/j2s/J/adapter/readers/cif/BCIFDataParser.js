Clazz.declarePackage("J.adapter.readers.cif");
Clazz.load(["JU.CifDataParser"], "J.adapter.readers.cif.BCIFDataParser", ["java.io.BufferedInputStream", "$.FileInputStream", "JU.BinaryDocument", "$.MessagePackReader", "J.adapter.readers.cif.BCIFDecoder"], function(){
var c$ = Clazz.decorateAsClass(function(){
this.rdr = null;
this.categoryName = null;
this.rowCount = 0;
this.columnMaps = null;
this.columnDcoders = null;
this.rowPt = 0;
this.ifield = 0;
this.dfield = 0;
this.fieldStr = null;
this.fieldIsValid = false;
this.cifMap = null;
this.header = null;
Clazz.instantialize(this, arguments);}, J.adapter.readers.cif, "BCIFDataParser", JU.CifDataParser);
Clazz.makeConstructor(c$, 
function(bcifReader, debugging){
Clazz.superConstructor (this, J.adapter.readers.cif.BCIFDataParser, []);
this.rdr = bcifReader;
this.debugging = debugging;
}, "J.adapter.readers.cif.BCIFReader,~B");
Clazz.defineMethod(c$, "getDecoder", 
function(key, col, rowCount, catName){
var sb = null;
var d =  new J.adapter.readers.cif.BCIFDecoder(key, col, sb).setRowCount(rowCount, catName).finalizeDecoding(null);
return (d == null || d.dataType == 0 ? null : d);
}, "~S,java.util.Map,~N,~S");
Clazz.defineMethod(c$, "initializeCategory", 
function(catName, rowCount, columns){
this.columnMaps = columns;
this.rowCount = rowCount;
this.categoryName = catName;
var n = this.columnCount = columns.length;
if (this.columnNames == null || this.columnNames.length < n) this.columnNames =  new Array(n);
for (var i = 0; i < n; i++) {
this.columnNames[i] = (catName + "_" + (columns[i]).get("name")).toLowerCase();
}
}, "~S,~N,~A");
Clazz.overrideMethod(c$, "parseDataBlockParameters", 
function(fieldNames, key, data, key2col, col2key){
this.haveData = false;
this.rowPt = -1;
if (!JU.CifDataParser.htFields.containsKey(fieldNames[0])) {
for (var i = fieldNames.length; --i >= 0; ) JU.CifDataParser.htFields.put(fieldNames[i], Integer.$valueOf(i));

}this.columnDcoders =  new Array(this.columnCount);
for (var pt = 0; pt < this.columnCount; pt++) {
var s = this.columnNames[pt];
var iField = JU.CifDataParser.htFields.get(s);
var keyIndex = col2key[pt] = (iField == null ? -1 : iField.intValue());
var d = (keyIndex == -1 ? null : this.getDecoder(null, this.getDataColumn(pt), this.rowCount, this.categoryName));
if (d == null) {
System.out.println("BCIFDataParser skipping " + s);
} else {
this.columnDcoders[pt] = d;
key2col[keyIndex] = pt;
this.haveData = true;
}}
}, "~A,~S,~S,~A,~A");
Clazz.defineMethod(c$, "decodeAndGetData", 
function(icol){
this.columnDcoders =  Clazz.newArray(-1, [this.getDecoder(null, this.getDataColumn(icol), this.rowCount, this.categoryName)]);
this.getColumnData(0);
this.columnDcoders = null;
}, "~N");
Clazz.overrideMethod(c$, "getData", 
function(){
this.rowPt++;
var done = this.rowPt >= this.rowCount;
if (done) {
for (var i = this.columnDcoders.length; --i >= 0; ) this.columnDcoders[i] = null;

}return !done;
});
Clazz.overrideMethod(c$, "getColumnData", 
function(colPt){
this.rdr.key = this.getColumnName(colPt);
this.ifield = -2147483648;
this.dfield = NaN;
if (this.columnDcoders[colPt] == null) {
this.fieldIsValid = false;
return this.fieldStr = this.nullString;
}switch (this.columnDcoders[colPt].dataType) {
case 1:
this.ifield = this.columnDcoders[colPt].getIntValue(this.rowPt);
this.fieldStr = (this.ifield == -2147483648 ? this.nullString : "_");
break;
case 2:
this.dfield = this.columnDcoders[colPt].getFixedPtValue(this.rowPt);
this.fieldStr = (Double.isNaN(this.dfield) ? this.nullString : "_");
break;
case 3:
this.fieldStr = this.columnDcoders[colPt].getStringValue(this.rowPt);
break;
}
this.fieldIsValid = (this.fieldStr !== this.nullString);
return this.fieldStr;
}, "~N");
Clazz.defineMethod(c$, "isFieldValid", 
function(){
return this.fieldIsValid;
});
Clazz.defineMethod(c$, "getDataColumn", 
function(icol){
return (icol >= 0 && icol < this.columnCount ? this.columnMaps[icol] : null);
}, "~N");
Clazz.overrideMethod(c$, "getAllCifData", 
function(){
return null;
});
c$.main = Clazz.defineMethod(c$, "main", 
function(args){
try {
var testFile = (args.length == 0 ? "c:/temp/1cbs.bcif" : args[0]);
var binaryDoc =  new JU.BinaryDocument();
var bis =  new java.io.BufferedInputStream( new java.io.FileInputStream(testFile));
binaryDoc.setStream(bis, true);
var msgMap;
msgMap = ( new JU.MessagePackReader(binaryDoc, false)).readMap();
var parser =  new J.adapter.readers.cif.BCIFDataParser(null, true);
parser.debugConstructCifMap(msgMap);
binaryDoc.close();
System.out.println("OK - DONE");
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
e.printStackTrace();
} else {
throw e;
}
}
}, "~A");
});
;//5.0.1-v2 Mon Feb 05 10:17:20 CST 2024
