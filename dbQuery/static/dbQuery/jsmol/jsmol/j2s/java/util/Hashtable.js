Clazz.load(["java.util.AbstractCollection", "$.AbstractSet", "$.Dictionary", "$.Enumeration", "$.Iterator", "$.Map"], "java.util.Hashtable", ["java.util.Collections"], function(){
var c$ = Clazz.decorateAsClass(function(){
this.table = null;
this.count = 0;
this.threshold = 0;
this.loadFactor = 0;
this.modCount = 0;
this.$keySet = null;
this.$entrySet = null;
this.$values = null;
this.秘m = null;
this.秘allowJS = false;
Clazz.instantialize(this, arguments);}, java.util, "Hashtable", java.util.Dictionary, [java.util.Map, Cloneable]);
Clazz.makeConstructor(c$, 
function(){
var map = null;
var capacity = 11;
var loadFactor = 0.75;
{
capacity = arguments[0];
loadFactor = arguments[1];
if (typeof capacity == "object") {
map = capacity;
capacity = Math.max(2*t.size(), 11);
this.秘allowJS = map.秘allowJS;
} else {
this.秘allowJS = true;
}
capacity = (capacity || 11);
loadFactor = (loadFactor || 0.75);
}if (capacity < 0) throw  new IllegalArgumentException("Illegal Capacity: " + capacity);
if (loadFactor <= 0 || Float.isNaN(loadFactor)) throw  new IllegalArgumentException("Illegal Load: " + loadFactor);
if (capacity == 0) capacity = 1;
this.loadFactor = loadFactor;
this.table =  new Array(capacity);
this.threshold = Clazz.floatToInt(Math.min(capacity * loadFactor, 2147483640));
this.秘setJS();
if (map != null) this.putAll(map);
});
Clazz.overrideMethod(c$, "size", 
function(){
var c = this.count;
{
c = this.秘m && this.秘m.size || c;
}return c;
});
Clazz.overrideMethod(c$, "isEmpty", 
function(){
return this.size() == 0;
});
Clazz.overrideMethod(c$, "keys", 
function(){
return this.getEnumeration(0);
});
Clazz.overrideMethod(c$, "elements", 
function(){
return this.getEnumeration(1);
});
Clazz.defineMethod(c$, "contains", 
function(value){
if (value == null) {
throw  new NullPointerException();
}if (this.size() == 0) return false;
if (java.util.Hashtable.秘isSimple(this)) {
var m = this.秘m;
{
var iter = m.values();
for (var n = iter.next(); !n.done; n = iter.next()) {
if (n.value == value || n.value.equals(value)) {
return true;
}
}
}} else {
var tab = this.table;
for (var i = tab.length; i-- > 0; ) {
for (var e = tab[i]; e != null; e = e.next_) {
if (e.value.equals(value)) {
return true;
}}
}
}return false;
}, "~O");
Clazz.overrideMethod(c$, "containsValue", 
function(value){
return this.contains(value);
}, "~O");
Clazz.overrideMethod(c$, "containsKey", 
function(key){
switch (java.util.Hashtable.秘hasKey(this, key)) {
case 0:
break;
case 1:
java.util.Hashtable.秘ensureJavaMap(this);
break;
case 2:
return false;
case 3:
return true;
}
var tab = this.table;
var hash = key.hashCode();
var index = (hash & 0x7FFFFFFF) % tab.length;
for (var e = tab[index]; e != null; e = e.next_) {
if ((e.hash == hash) && e.key.equals(key)) {
return true;
}}
return false;
}, "~O");
Clazz.overrideMethod(c$, "get", 
function(key){
if (key == null) return null;
switch (java.util.Hashtable.秘hasKey(this, key)) {
case 0:
break;
case 1:
java.util.Hashtable.秘ensureJavaMap(this);
break;
case 2:
return null;
case 3:
var v = null;
{
v = this.秘m.get(key);
}return v;
}
var tab = this.table;
var hash = key.hashCode();
var index = (hash & 0x7FFFFFFF) % tab.length;
for (var e = tab[index]; e != null; e = e.next_) {
if ((e.hash == hash) && e.key.equals(key)) {
return e.value;
}}
return null;
}, "~O");
Clazz.defineMethod(c$, "rehash", 
function(){
var oldCapacity = this.table.length;
var oldMap = this.table;
var newCapacity = (oldCapacity << 1) + 1;
if (newCapacity - 2147483639 > 0) {
if (oldCapacity == 2147483639) return;
newCapacity = 2147483639;
}var newMap =  new Array(newCapacity);
this.modCount++;
this.threshold = Clazz.floatToInt(Math.min(newCapacity * this.loadFactor, 2147483640));
this.table = newMap;
for (var i = oldCapacity; i-- > 0; ) {
for (var old = oldMap[i]; old != null; ) {
var e = old;
old = old.next_;
var index = (e.hash & 0x7FFFFFFF) % newCapacity;
e.next_ = newMap[index];
newMap[index] = e;
}
}
});
Clazz.defineMethod(c$, "addEntry", 
function(hash, key, value, index){
this.modCount++;
var tab = this.table;
if (this.count >= this.threshold) {
this.rehash();
tab = this.table;
hash = key.hashCode();
index = (hash & 0x7FFFFFFF) % tab.length;
}var e = tab[index];
tab[index] =  new java.util.Hashtable.Entry(hash, key, value, e);
this.count++;
}, "~N,~O,~O,~N");
Clazz.overrideMethod(c$, "put", 
function(key, value){
if (value == null) {
throw  new NullPointerException();
}switch (java.util.Hashtable.秘hasKey(this, key)) {
case 0:
break;
case 1:
java.util.Hashtable.秘ensureJavaMap(this);
break;
case 2:
{
this.秘m.set(key, value);
}++this.modCount;
return null;
case 3:
var v0 = null;
{
v0 = this.秘m.get(key);
this.秘m.set(key, value);
}++this.modCount;
return v0;
}
var tab = this.table;
var hash = key.hashCode();
var index = (hash & 0x7FFFFFFF) % tab.length;
var entry = tab[index];
for (; entry != null; entry = entry.next_) {
if ((entry.hash == hash) && entry.key.equals(key)) {
var old = entry.value;
entry.value = value;
return old;
}}
this.addEntry(hash, key, value, index);
return null;
}, "~O,~O");
Clazz.overrideMethod(c$, "remove", 
function(key){
if (key == null) throw  new NullPointerException("Hashtable key may not be null");
switch (java.util.Hashtable.秘hasKey(this, key)) {
case 0:
break;
case 1:
java.util.Hashtable.秘ensureJavaMap(this);
break;
case 2:
return null;
case 3:
var v0 = null;
{
v0 = this.秘m.get(key); this.秘m["delete"](key);
}++this.modCount;
return v0;
}
var tab = this.table;
var hash = key.hashCode();
var index = (hash & 0x7FFFFFFF) % tab.length;
var e = tab[index];
for (var prev = null; e != null; prev = e, e = e.next_) {
if ((e.hash == hash) && e.key.equals(key)) {
this.modCount++;
if (prev != null) {
prev.next_ = e.next_;
} else {
tab[index] = e.next_;
}this.count--;
var oldValue = e.value;
e.value = null;
return oldValue;
}}
return null;
}, "~O");
Clazz.overrideMethod(c$, "putAll", 
function(t){
var key = null;
var value = null;
if (java.util.Hashtable.秘isSimple(t)) {
var me = this;
{
t.秘m.forEach(function(value, key) { me.put(key, value); })
}return;
}for (var e, $e = t.entrySet().iterator (); $e.hasNext()&& ((e = $e.next ()) || true);) this.put(e.getKey(), e.getValue());

}, "java.util.Map");
Clazz.overrideMethod(c$, "clear", 
function(){
var tab = this.table;
this.modCount++;
if (java.util.Hashtable.秘isSimple(this)) {
{
this.秘m.clear();
}}this.秘setJS();
for (var index = tab.length; --index >= 0; ) tab[index] = null;

this.count = 0;
});
Clazz.defineMethod(c$, "clone", 
function(){
try {
var t = Clazz.superCall(this, java.util.Hashtable, "clone", []);
t.table =  new Array(this.table.length);
for (var i = this.table.length; i-- > 0; ) {
t.table[i] = (this.table[i] != null) ? this.table[i].clone() : null;
}
t.$keySet = null;
t.$entrySet = null;
t.$values = null;
t.modCount = 0;
if (java.util.Hashtable.秘isSimple(this)) {
t.秘setJS();
var me = this;
{
me.秘m.forEach(function(value, key) {
t.秘m.set(key, value); t.modCount++;
});
}} else {
t.秘m = null;
}return t;
} catch (e) {
if (Clazz.exceptionOf(e,"CloneNotSupportedException")){
throw  new InternalError(e);
} else {
throw e;
}
}
});
Clazz.defineMethod(c$, "toString", 
function(){
var max = this.size() - 1;
if (max == -1) return "{}";
var sb =  new StringBuilder();
var it = this.entrySet().iterator();
sb.append('{');
for (var i = 0; ; i++) {
var e = it.next();
var key = e.getKey();
var value = e.getValue();
sb.append(key === this ? "(this Map)" : key.toString());
sb.append('=');
sb.append(value === this ? "(this Map)" : value.toString());
if (i == max) return sb.append('}').toString();
sb.append(", ");
}
});
Clazz.overrideMethod(c$, "equals", 
function(o){
if (o === this) return true;
if (!(Clazz.instanceOf(o,"java.util.Map"))) return false;
var t = o;
if (t.size() != this.size()) return false;
try {
var i = this.entrySet().iterator();
while (i.hasNext()) {
var e = i.next();
var key = e.getKey();
var value = e.getValue();
if (value == null) {
if (!(t.get(key) == null && t.containsKey(key))) return false;
} else {
if (!value.equals(t.get(key))) return false;
}}
} catch (e$$) {
if (Clazz.exceptionOf(e$$,"ClassCastException")){
var unused = e$$;
{
return false;
}
} else if (Clazz.exceptionOf(e$$, NullPointerException)){
var unused = e$$;
{
return false;
}
} else {
throw e$$;
}
}
return true;
}, "~O");
Clazz.overrideMethod(c$, "hashCode", 
function(){
var h = 0;
if (this.count == 0 || this.loadFactor < 0) return h;
this.loadFactor = -this.loadFactor;
var tab = this.table;
for (var entry, $entry = 0, $$entry = tab; $entry < $$entry.length && ((entry = $$entry[$entry]) || true); $entry++) {
while (entry != null) {
h += entry.hashCode();
entry = entry.next_;
}
}
this.loadFactor = -this.loadFactor;
return h;
});
Clazz.defineMethod(c$, "getEnumeration", 
function(type){
if (this.size() == 0) {
return java.util.Collections.emptyEnumeration();
} else {
return  new java.util.Hashtable.Enumerator(this, type, false);
}}, "~N");
Clazz.defineMethod(c$, "getIterator", 
function(type){
if (this.size() == 0) {
return java.util.Collections.emptyIterator();
} else {
return  new java.util.Hashtable.Enumerator(this, type, true);
}}, "~N");
Clazz.overrideMethod(c$, "keySet", 
function(){
if (this.$keySet == null) this.$keySet =  new java.util.Hashtable.KeySet(this);
return this.$keySet;
});
Clazz.overrideMethod(c$, "entrySet", 
function(){
if (this.$entrySet == null) this.$entrySet =  new java.util.Hashtable.EntrySet(this);
return this.$entrySet;
});
Clazz.overrideMethod(c$, "values", 
function(){
if (this.$values == null) this.$values =  new java.util.Hashtable.ValueCollection(this);
return this.$values;
});
Clazz.defineMethod(c$, "秘setJS", 
function(){
if (this.秘allowJS && java.util.Hashtable.USE_SIMPLE) {
var m = null;
{
m = new Map();
}this.秘m = m;
} else {
this.秘m = null;
}});
c$.秘get = Clazz.defineMethod(c$, "秘get", 
function(map, key){
{
return map.秘m.get(key == null ? null : key + "")
}}, "~O,~O");
c$.秘set = Clazz.defineMethod(c$, "秘set", 
function(map, key, value){
}, "java.util.Map,~O,~O");
c$.秘hasKey = Clazz.defineMethod(c$, "秘hasKey", 
function(map, key){
{
return (!map.秘m ? 0 : key != null && typeof key != "string"
? 1 : map.秘m.has(key) ? 3 : 2);
}}, "java.util.Map,~O");
c$.秘isSimple = Clazz.defineMethod(c$, "秘isSimple", 
function(map){
{
return !!map.秘m;
}}, "java.util.Map");
c$.秘ensureJavaMap = Clazz.defineMethod(c$, "秘ensureJavaMap", 
function(map){
}, "java.util.Map");
/*if3*/;(function(){
var c$ = Clazz.decorateAsClass(function(){
this.hash = 0;
this.key = null;
this.value = null;
this.next_ = null;
Clazz.instantialize(this, arguments);}, java.util.Hashtable, "Entry", null, java.util.Map.Entry);
Clazz.makeConstructor(c$, 
function(hash, key, value, next){
this.hash = hash;
this.key = key;
this.value = value;
this.next_ = next;
}, "~N,~O,~O,java.util.Hashtable.Entry");
Clazz.overrideMethod(c$, "clone", 
function(){
return  new java.util.Hashtable.Entry(this.hash, this.key, this.value, (this.next_ == null ? null : this.next_.clone()));
});
Clazz.overrideMethod(c$, "getKey", 
function(){
return this.key;
});
Clazz.overrideMethod(c$, "getValue", 
function(){
return this.value;
});
Clazz.overrideMethod(c$, "setValue", 
function(value){
if (value == null) throw  new NullPointerException();
var oldValue = this.value;
this.value = value;
return oldValue;
}, "~O");
Clazz.overrideMethod(c$, "equals", 
function(o){
if (!(Clazz.instanceOf(o,"java.util.Map.Entry"))) return false;
var e = o;
return (this.key == null ? e.getKey() == null : this.key.equals(e.getKey())) && (this.value == null ? e.getValue() == null : this.value.equals(e.getValue()));
}, "~O");
Clazz.overrideMethod(c$, "hashCode", 
function(){
return this.hash ^ (this.value == null ? 0 : this.value.hashCode());
});
Clazz.defineMethod(c$, "toString", 
function(){
return this.key.toString() + "=" + this.value.toString();
});
/*eoif3*/})();
/*if3*/;(function(){
var c$ = Clazz.decorateAsClass(function(){
this.ht = null;
Clazz.instantialize(this, arguments);}, java.util.Hashtable, "KeySet", java.util.AbstractSet);
Clazz.makeConstructor(c$, 
function(ht){
Clazz.superConstructor (this, java.util.Hashtable.KeySet, []);
this.ht = ht;
}, "java.util.Hashtable");
Clazz.overrideMethod(c$, "iterator", 
function(){
return this.ht.getIterator(0);
});
Clazz.overrideMethod(c$, "size", 
function(){
return this.ht.size();
});
Clazz.overrideMethod(c$, "contains", 
function(o){
return this.ht.containsKey(o);
}, "~O");
Clazz.overrideMethod(c$, "remove", 
function(o){
return this.ht.remove(o) != null;
}, "~O");
Clazz.overrideMethod(c$, "clear", 
function(){
this.ht.clear();
});
/*eoif3*/})();
/*if3*/;(function(){
var c$ = Clazz.decorateAsClass(function(){
this.ht = null;
Clazz.instantialize(this, arguments);}, java.util.Hashtable, "EntrySet", java.util.AbstractSet);
Clazz.makeConstructor(c$, 
function(ht){
Clazz.superConstructor (this, java.util.Hashtable.EntrySet, []);
this.ht = ht;
}, "java.util.Hashtable");
Clazz.overrideMethod(c$, "iterator", 
function(){
return this.ht.getIterator(2);
});
Clazz.overrideMethod(c$, "contains", 
function(o){
if (!(Clazz.instanceOf(o,"java.util.Map.Entry"))) return false;
var entry = o;
var key = entry.getKey();
switch (java.util.Hashtable.秘hasKey(this.ht, key)) {
case 0:
break;
case 1:
java.util.Hashtable.秘ensureJavaMap(this.ht);
break;
case 3:
var value = entry.getValue();
var v = this.ht.get(key);
return (value === v || value != null && value.equals(key));
case 2:
return false;
}
var tab = this.ht.table;
var hash = key.hashCode();
var index = (hash & 0x7FFFFFFF) % tab.length;
for (var e = tab[index]; e != null; e = e.next_) if (e.hash == hash && e.equals(entry)) return true;

return false;
}, "~O");
Clazz.overrideMethod(c$, "remove", 
function(o){
if (!(Clazz.instanceOf(o,"java.util.Map.Entry"))) return false;
var entry = o;
var key = entry.getKey();
switch (java.util.Hashtable.秘hasKey(this.ht, key)) {
case 0:
break;
case 1:
java.util.Hashtable.秘ensureJavaMap(this.ht);
break;
case 3:
var value = entry.getValue();
if (value == null) return false;
var v = this.ht.get(key);
if (v === value || v.equals(value)) {
this.ht.remove(key);
return true;
}return false;
case 2:
return false;
}
var tab = this.ht.table;
var hash = key.hashCode();
var index = (hash & 0x7FFFFFFF) % tab.length;
var e = tab[index];
for (var prev = null; e != null; prev = e, e = e.next_) {
if (e.hash == hash && e.equals(entry)) {
this.ht.modCount++;
if (prev != null) prev.next_ = e.next_;
 else tab[index] = e.next_;
this.ht.count--;
e.value = null;
return true;
}}
return false;
}, "~O");
Clazz.overrideMethod(c$, "size", 
function(){
return this.ht.size();
});
Clazz.overrideMethod(c$, "clear", 
function(){
this.ht.clear();
});
/*eoif3*/})();
/*if3*/;(function(){
var c$ = Clazz.decorateAsClass(function(){
this.ht = null;
Clazz.instantialize(this, arguments);}, java.util.Hashtable, "ValueCollection", java.util.AbstractCollection);
Clazz.makeConstructor(c$, 
function(ht){
Clazz.superConstructor (this, java.util.Hashtable.ValueCollection, []);
this.ht = ht;
}, "java.util.Hashtable");
Clazz.overrideMethod(c$, "iterator", 
function(){
return this.ht.getIterator(1);
});
Clazz.overrideMethod(c$, "size", 
function(){
return this.ht.size();
});
Clazz.overrideMethod(c$, "contains", 
function(o){
return this.ht.containsValue(o);
}, "~O");
Clazz.overrideMethod(c$, "clear", 
function(){
this.ht.clear();
});
/*eoif3*/})();
/*if3*/;(function(){
var c$ = Clazz.decorateAsClass(function(){
this.table = null;
this.index = 0;
this.next_ = null;
this.current = null;
this.type = 0;
this.jsMapIterator = null;
this.ht = null;
this.isIterator = false;
this.expectedModCount = 0;
Clazz.instantialize(this, arguments);}, java.util.Hashtable, "Enumerator", null, [java.util.Enumeration, java.util.Iterator]);
Clazz.makeConstructor(c$, 
function(ht, type, iterator){
this.ht = ht;
this.table = ht.table;
this.index = ht.table.length;
this.type = type;
this.isIterator = iterator;
this.expectedModCount = ht.modCount;
if (java.util.Hashtable.秘isSimple(ht)) {
var m = ht.秘m;
{
this.jsMapIterator = m.entries();
this.next_ = this.jsMapIterator.next();
}}}, "java.util.Hashtable,~N,~B");
Clazz.overrideMethod(c$, "hasMoreElements", 
function(){
if (java.util.Hashtable.秘isSimple(this.ht)) {
var b = false;
{
b = this.next_ && !this.next_.done
}return b;
} else {
var e = this.next_;
var i = this.index;
var t = this.table;
while (e == null && i > 0) {
e = t[--i];
}
this.next_ = e;
this.index = i;
return e != null;
}});
Clazz.overrideMethod(c$, "nextElement", 
function(){
var node = this.next_;
if (java.util.Hashtable.秘isSimple(this.ht)) {
var t = this.type;
this.current = node;
var n = null;
{
n = this.jsMapIterator.next() || null;
}this.next_ = n;
if (node != null) {
var k = null;
var v = null;
var done = false;
{
done = node.done;
if (!done) {
if (t < 2) return node.value[t];
k = node.value[0];
v = node.value[1];
}
}if (!done) {
return ((Clazz.isClassDefined("java.util.Hashtable$Enumerator$1") ? 0 : java.util.Hashtable.Enumerator.$Hashtable$Enumerator$1$ ()), Clazz.innerTypeInstance(java.util.Hashtable$Enumerator$1, this, null, 0, k, v, null));
}}} else {
var i = this.index;
var t = this.table;
while (node == null && i > 0) {
node = t[--i];
}
this.next_ = node;
this.index = i;
if (node != null) {
var e = this.current = this.next_;
this.next_ = e.next_;
return this.type == 0 ? e.key : (this.type == 1 ? e.value : e);
}}throw  new java.util.NoSuchElementException("Hashtable Enumerator");
});
Clazz.overrideMethod(c$, "hasNext", 
function(){
return this.hasMoreElements();
});
Clazz.overrideMethod(c$, "next", 
function(){
if (this.ht.modCount != this.expectedModCount) throw  new java.util.ConcurrentModificationException();
return this.nextElement();
});
Clazz.overrideMethod(c$, "remove", 
function(){
if (!this.isIterator) throw  new UnsupportedOperationException();
var p = this.current;
if (p == null) throw  new IllegalStateException("Hashtable Enumerator");
if (this.ht.modCount != this.expectedModCount) throw  new java.util.ConcurrentModificationException();
if (java.util.Hashtable.秘isSimple(this.ht)) {
var key = null;
{
key = p.value[0];
}this.ht.remove(key);
this.expectedModCount++;
} else {
{
var tab = this.ht.table;
var index = (this.current.hash & 0x7FFFFFFF) % tab.length;
var e = tab[index];
for (var prev = null; e != null; prev = e, e = e.next_) {
if (e === this.current) {
this.ht.modCount++;
this.expectedModCount++;
if (prev == null) tab[index] = e.next_;
 else prev.next_ = e.next_;
this.ht.count--;
this.current = null;
return;
}}
throw  new java.util.ConcurrentModificationException();
}}});
c$.$Hashtable$Enumerator$1$=function(){
/*if5*/;(function(){
var c$ = Clazz.declareAnonymous(java.util, "Hashtable$Enumerator$1", java.util.Hashtable.Entry);
Clazz.overrideMethod(c$, "setValue", 
function(value){
var m = this.b$["java.util.Hashtable.Enumerator"].ht.modCount;
var v = this.b$["java.util.Hashtable.Enumerator"].ht.put(this.getKey(), value);
this.b$["java.util.Hashtable.Enumerator"].ht.modCount = m;
return v;
}, "~O");
/*eoif5*/})();
};
/*eoif3*/})();
c$.USE_SIMPLE = true;
});
;//5.0.1-v2 Mon Jan 22 13:18:57 CST 2024
