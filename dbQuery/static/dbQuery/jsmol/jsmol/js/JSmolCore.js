// JSmolCore.js -- Jmol core capability 

// allows Jmol applets to be created on a page with more flexibility and extendability
// provides an object-oriented interface for JSpecView and syncing of Jmol/JSpecView

// see JSmolApi.js for public user-interface. All these are private functions

// BH 2023.01.19 Jmol._allowKeyboardFocus = true; can be set to false if page jumping is a problem.
// BH 8/25/2022 fixes getFileData response "" from empty file returning "OK" instead of ""
// BH 8/15/2022 adds .lut for binary
// BH 6/23/2022 implements Jmol._lastAppletID via setMouseOwner
// BH 5/12/2022 adds setting file type option for drag drop
// BH 4/30/2019 fixes write xyz "https://...."
// BH 7/6/2017 2:22:07 AM adds BZ2 as binary
// BH 4/13/2017 11:23:05 PM adds "binary pmesh" .pmb extension
// BH 1/14/2017 6:28:07 AM adds &debugCore
// BH 10/20/2016 10:00:43 AM JmolTracker.php
// BH 9/19/2016 8:22:48 AM drag-drop broken for https (imageDrop.htm)
// BH 9/18/2016 btoa() does not work with UTF-8 data (set language es;write menu t.mnu)
// BH 8/26/2016 11:29:48 AM RCSB ligand 
// BH 8/26/2016 11:29:48 AM generic fixProtocol for .gov/ to https
// BH 8/26/2016 6:56:31 AM chemapps.stolaf.edu exclusively https
// BH 8/25/2016 9:47:26 PM bug fix: NCI/CADD now requires "get3d=true" not "get3d=True"
// BH 7/31/2016 6:42:06 AM changes mouse wheel from -1 to 507
// BH 6/27/2016 1:16:57 AM adds Jmol.playAudio(fname)
// BH 4/26/2016 4:16:07 PM adds Jmol.loadFileFromDialog(applet)
// BH 4/21/2016 9:25:39 AM adds [URL] button to file load option
// BH 4/20/2016 2:44:50 PM fixes async load problem with Safari
// BH 4/18/2016 10:25:08 PM adds preliminary =xxxx.mmtf reader
// BH 4/13/2016 9:12:31 PM  url.indexOf("http://www.rcsb.org/pdb/files/") == 0 && url.indexOf("/ligand/") < 0 ? 
// BH 4/11/2016 5:34:16 PM adds direct conversion to http://files.rcsb.org/view from http://www.rcsb.org/pdb/files/1xpb.pdb
// BH 4/3/2016 9:10:31 PM adding materialsproject.org for AJAX.
// BH 3/23/2016 1:21:39 PM adding http://files.rcsb.org/view/%FILE.pdb as default RCSB site for "="

// BH 2/29/2016 3:59:55 PM broken cursor_wait image path when Info.j2sPath is not "j2s"
// BH 2/19/2016 10:32:18 AM typo fixed for makeLiveImage
// BH 2/14/2016 12:31:02 PM fixed local reader not disappearing after script call
// BH 2/14/2016 12:30:41 PM Info.appletLoadingImage: "j2s/img/JSmol_spinner.gif", 
   // can be set to "none" or some other image; see Jmol._hideLoadingSpinner(applet)
   // implemented only for JSmolApplet, not others
// BH 2/14/2016 12:27:09 PM Jmol.setCursor 
// BH 2/14/2016 6:48:33 AM _setCursor() and cursor_wait   http://ajaxload.info/
// BH 1/15/2016 4:23:14 PM adding Info.makeLiveImage
// BH 12/30/2015 8:18:42 PM adding AMS call to database list; allowing for ?ALLOWSORIGIN? to override settings here
// BH 12/17/2015 4:43:05 PM adding Jmol._requestRepaint to allow for MSIE9 not having requestAnimationFrame
// BH 12/16/2015 3:01:06 PM adding $.ajaxSetup({ mimeType: "text/plain" });
// BH 12/14/2015 6:42:03 PM adding check for MS Edge browser, which does not support dataURI
// BH 12/2/2015 1:18:15 PM adding .dcd as binary file type
// BH 12/1/2015 10:05:55 AM loading identical HTML5 page after Java page causes bad NPObject error 
// BH 10/26/2015 12:47:16 PM adding two rcsb sites for direct access
// BH 10/23/2015 9:20:39 PM minor coding adjustment
// BH 10/13/2015 9:32:08 PM adding Jmol.__$ as jquery object used 
// BH 15/09/2015 18:06:39 fixing mouse check for swingjs-ui since SVG element className is not a string 
// BH 8/12/2015 11:43:52 PM adding isHttps2Http forcing call to server proxy
// BH 8/9/2015 6:33:33 PM correcting bug in load ASYNC for x-domain access
// BH 7/7/2015 1:42:31 PM Jmol._persistentMenu
// BH 6/29/2015 10:14:47 AM adds Jmol.$getSize(obj)
// BH 5/30/2015 9:33:12 AM adds class swingjs-ui to ignore 
// BH 5/9/2015 3:38:52 PM adds data-ignoreMouse attribute for JTextField
// BH 3/30/2015 9:46:53 PM adds JSAppletPanel for ready callback
// BH 12/6/2014 3:32:54 PM Jmol.setAppletCss() broken
// BH 9/13/2014 2:15:51 PM embedded JSME loads from SEARCH when Jmol should 
// BH 8/14/2014 2:52:38 PM drag-drop cache should not be cleared if SPT file is dropped
// BH 8/5/2014 6:39:54 AM unnecessary messages about binary for PDB finally removed
// BH 8/4/2014 5:30:00 AM automatically switch to no document after page loading
// BH 8/2/2014 5:22:40 PM drag-drop broken in JSmol/HTML5 
// BH 7/23/2014 5:34:08 PM setting a parameter such as readyFunction to null stops file loading
// BH 7/3/2014 12:30:28 AM lost drag-drop of models
// BH 7/2/2014 4:47:55 AM adding pdbe.org to direct database calls
// BH 5/30/2014 7:20:07 AM better dragging for console and menu
// BH 4/27/2014 6:31:52 PM allows _USE=SIGNED HTML5 as well as _USE=JAVA HTML5
// BH 3/8/2014 5:50:51 PM adds support for dataURI download in FF and Chrome
// BH 3/8/2014 8:43:10 AM moves PubChem access to https
// BH 3/4/2014 8:40:15 PM adds Jmol.Cache for JSV/Jmol sharing files
// BH 2/10/2014 10:07:14 AM added Info.z and Info.zIndexBase
// BH 2/9/2014 9:56:06 PM updated JSmolCore.js with option to extend Viewer with code PRIOR to loading Viewer classes
// BH 2/6/2014 8:46:25 AM disabled Jmol._tracker for localhost and 127.x.x.x 
// BH 1/29/2014 8:02:23 AM Jmol.View and Info.viewSet
// BH 1/21/2014 12:06:59 PM adding Jmol.Info.cacheFiles (applet, true/false) and applet._cacheFiles and Jmol._fileCache
// BH 1/13/2014 2:12:38 PM adding "http://www.nmrdb.org/tools/jmol/predict.php":"%URL", to _DirectDatabaseCalls
// BH 12/21/2013 6:38:35 PM applet sync broken
// BH 12/6/2013 6:18:32 PM cover.htm and coverImage fix
// BH 12/4/2013 7:44:26 PM fix for JME independent search box
// BH 12/3/2013 6:30:08 AM fix for ready function returning Boolean instead of boolean in HTML5 version
// BH 11/30/2013 10:31:37 AM added type:"GET" for jQuery.ajax() requests instead of using defaults
// BH 11/30/2013 10:31:37 AM added cache:true for jQuery.ajax() requests; can override with cache:"NO", not cache:false
// BH 11/28/2013 11:09:27 AM added Jmol._alertNoBinary:true
// BH 11/26/2013 8:19:55 PM fix !xxxx search commmand entry and stop MSIE from duplicating command
// BH 11/25/2013 7:38:31 AM adds Jmol._tracker: option for GoogleAnalytics tracking
// BH 11/25/2013 7:39:03 AM adds URL options _J2S=  _JAR=  _USE=
// BH 11/23/2013 10:51:37 PM  adds JNLP support for local applet
// BH 11/2/2013 12:05:11 PM JSmolJSME fixes; https access fixed
// BH 10/31/2013 7:50:06 PM Jmol.Dialog as SwingController; Jmol._mouseOwner added
// BH 10/19/2013 7:05:04 AM adding Jmol._ajaxCall for Error Contacting Server; database POST method enabled
// BH 10/17/2013 1:40:51 PM  adding javajs/swing and Jmol.Dialog
// BH 9/30/2013 6:42:24 PM: pdb.gz switch  pdb should only be for www.rcsb.org
// BH 9/17/2013 10:17:51 AM: asynchronous file reading and saving
// BH 8/16/2013 12:02:20 PM: JSmoljQueryExt.js pulled out
// BH 8/16/2013 12:02:20 PM: Jmol._touching used properly

// BH 3/22/2013 5:53:02 PM: Adds noscript option, JSmol.min.core.js
// BH 1/17/2013 5:20:44 PM: Fixed problem with console not getting initial position if no first click
// 1/13/2013 BH: Fixed MSIE not-reading-local-files problem.
// 11/28/2012 BH: Fixed MacOS Safari binary ArrayBuffer problem
// 11/21/2012 BH: restructuring of files as JS... instead of J...
// 11/20/2012 BH: MSIE9 cannot do a synchronous file load cross-domain. See Jmol._getFileData
// 11/4/2012 BH: RCSB REST format change "<structureId>" to "<dimStructure.structureId>"
// 9/13/2012 BH: JmolCore.js changes for JSmol doAjax() method -- _3ata()
// 6/12/2012 BH: JmolApi.js: adds Jmol.setInfo(applet, info, isShown) -- third parameter optional 
// 6/12/2012 BH: JmolApi.js: adds Jmol.getInfo(applet) 
// 6/12/2012 BH: JmolApplet.js: Fixes for MSIE 8
// 6/5/2012  BH: fixes problem with Jmol "javascript" command not working and getPropertyAsArray not working
// 6/4/2012  BH: corrects problem with MSIE requiring mouse-hover to activate applet
// 5/31/2012 BH: added JSpecView interface and api -- see JmolJSV.js
//               also changed "jmolJarPath" to just "jarPath"
//               jmolJarFile->jarFile, jmolIsSigned->isSigned, jmolReadyFunction->readyFunction
//               also corrects a double-loading issue
// 5/14/2012 BH: added AJAX queue for ChemDoodle option with multiple canvases 
// 8/12/2012 BH: adds support for MSIE xdr cross-domain request (jQuery.iecors.js)

	// BH 4/25 -- added text option. setAppletCss(null, "style=\"xxxx\"")
	// note that since you must add the style keyword, this can be used to add any attribute to these tags, not just css. 

// required/optional libraries (preferably in the following order):

//    jquery/jquery.js     -- at least jQuery.1.9
//    js/JSmoljQueryext.js -- required for binary file transfer; otherwise standard jQuery should be OK
//    js/JSmolCore.js      -- required
//    js/j2sjmol.js        -- required
//    js/JSmol.js          -- required
//    js/JSmolApplet.js    -- required; internal functions for _Applet and _Image; must be after JSmolCore
//    js/JSmolControls.js  -- optional; internal functions for buttons, links, menus, etc.; must be after JSmolCore
//    js/JSmolConsole.js   -- optional; for the pop-up console
//    js/JSmolApi.js       -- required; all user functions; must be after JSmolCore
//    js/JSmolTHREE.js     -- optional; WebGL library required for JSmolGLmol.js
//    js/JSmolGLmol.js     -- optional; WebGL version of JSmol.
//    js/JSmolJME.js       -- optional; JSME (2D editor)
//    jsme/jsme/jsme.nocache.js   --  required for JSME 
//    js/JSmolMenu.js      -- optional; required for menuing in JSV
//    js/JSmolJSV.js       -- optional; for creating and interacting with a JSpecView applet 

// most of these will be loaded automatically, and for most installations, all you need is JSmol.min.js


// Allows Jmol-like objects to be displayed on Java-challenged (iPad/iPhone)
// or applet-challenged (Android/iPhone) platforms, with automatic switching to 

// For your installation, you should consider putting JmolData.jar and jsmol.php 
// on your own server. Nothing more than these two files is needed on the server, and this 
// allows more options for MSIE and Chrome when working with cross-domain files (such as RCSB or pubChem) 

// The NCI and RCSB databases are accessed via direct AJAX if available (xhr2/xdr).


if(typeof(jQuery)=="undefined") alert ("Note -- JSmoljQuery is required for JSmol, but it's not defined.")

// An example of how to extend Jmol with code PRIOR to JSmolCore.js or JSmol.min.js follows:
//
// 
//	Jmol = {
//  	z:3000,
//		extend: function(what, obj) {if (what == "viewer") { obj._testing = true } }
//	}

self.Jmol || (Jmol = {});

if (!Jmol._version)
Jmol = (function(document) {
	var z=Jmol.z || 9000;
	var getZOrders = function(z) {
		return {
			rear:z++,
			header:z++,
			main:z++,
			image:z++,
			front:z++,
			fileOpener:z++,
			coverImage:z++,
			dialog:z++, // could be several of these, JSV only
			menu:z+90000, // way front
			console:z+91000, // even more front
      consoleImage:z+91001, // bit more front; increments
			monitorZIndex:z+99999 // way way front
		}
	};
	var j = {
		_version: "$Date: 2022-06-24 05:54:49 -0500 (Fri, 24 Jun 2022) $", // svn.keywords:lastUpdated
		_alertNoBinary: true,
		// this url is used to Google Analytics tracking of Jmol use. You may remove it or modify it if you wish. 
		_allowedJmolSize: [25, 2048, 300],   // min, max, default (pixels)
		/*  By setting the Jmol.allowedJmolSize[] variable in the webpage
				before calling Jmol.getApplet(), limits for applet size can be overriden.
				2048 standard for GeoWall (http://geowall.geo.lsa.umich.edu/home.html)
		*/
		_allowKeyboardFocus: true, // set false on page if mouse moves and clicks cause too much scrolling of window
		_appletCssClass: "",
		_appletCssText: "",
		_fileCache: null, // enabled by Jmol.setFileCaching(applet, true/false)
		_jarFile: null,  // can be set in URL using _JAR=
		_j2sPath: null,  // can be set in URL using _J2S=
		_use: null,      // can be set in URL using _USE=
		_j2sLoadMonitorOpacity: 90, // initial opacity for j2s load monitor message
		_applets: {},
		_asynchronous: true,
		_ajaxQueue: [],
    _persistentMenu: false,
		_getZOrders: getZOrders,
		_z:getZOrders(z),
		_debugCode: true,  // set false in process of minimization
    _debugCore: false, // set true using URL &debugCore
		db: {
			_databasePrefixes: "$=:",
			_fileLoadScript: ";if (_loadScript = '' && defaultLoadScript == '' && _filetype == 'Pdb') { select protein or nucleic;cartoons Only;color structure; select * };",
			_nciLoadScript: ";n = ({molecule=1}.length < {molecule=2}.length ? 2 : 1); select molecule=n;display selected;center selected;",
			_pubChemLoadScript: "",
			_DirectDatabaseCalls:{
				// these sites are known to implement access-control-allow-origin *
        // null here means no conversion necessary 
				"cactus.nci.nih.gov": null,
        ".x3dna.org": null,
        "rruff.geo.arizona.edu": null, 
        ".rcsb.org": null, 
				"ftp.wwpdb.org": null,
				"pdbe.org": null, 
				"materialsproject.org": null, 
				".ebi.ac.uk": null, 
				"pubchem.ncbi.nlm.nih.gov":null,
				"www.nmrdb.org/tools/jmol/predict.php":null,
				"$": "https://cactus.nci.nih.gov/chemical/structure/%FILENCI/file?format=sdf&get3d=True",
				"$$": "https://cactus.nci.nih.gov/chemical/structure/%FILENCI/file?format=sdf",
				"=": "https://files.rcsb.org/download/%FILE.pdb",
				"*": "https://www.ebi.ac.uk/pdbe/entry-files/download/%FILE.cif",
				"==": "https://files.rcsb.org/ligands/download/%FILE.cif",
				":": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/%FILE/SDF?record_type=3d"
			},
			_restQueryUrl: "http://www.rcsb.org/pdb/rest/search",
			_restQueryXml: "<orgPdbQuery><queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType><description>Text Search</description><keywords>QUERY</keywords></orgPdbQuery>",
			_restReportUrl: "http://www.pdb.org/pdb/rest/customReport?pdbids=IDLIST&customReportColumns=structureId,structureTitle"
		},
		_debugAlert: false,
		_document: document,
		_isXHTML: false,
		_lastAppletID: null,
		_mousePageX: null,
		_mouseOwner: null,
		_serverUrl: "https://your.server.here/jsmol.php",
		_syncId: ("" + Math.random()).substring(3),
		_touching: false,
		_XhtmlElement: null,
		_XhtmlAppendChild: false
	}
	
	var ref = document.location.href.toLowerCase();
	j._debugCore = (ref.indexOf("j2sdebugcore") >= 0);
	// j._debugCode is set by the JmolDebugOff.js script fragment
	j._httpProto = (ref.indexOf("https") == 0 ? "https://" : "http://"); 
	j._isFile = (ref.indexOf("file:") == 0);
	if (j._isFile) // ensure no attempt to read XML in local request:
	  $.ajaxSetup({ mimeType: "text/plain" });
	j._ajaxTestSite = j._httpProto + "google.com";
	var isLocal = (j._isFile || ref.indexOf("http://localhost") == 0 || ref.indexOf("http://127.") == 0);
	j._tracker = (!isLocal && 'https://chemapps.stolaf.edu/jmol/JmolTracker.php?id=UA-45940799-1');
	
	j._isChrome = (navigator.userAgent.toLowerCase().indexOf("chrome") >= 0);
	j._isSafari = (!j._isChrome && navigator.userAgent.toLowerCase().indexOf("safari") >= 0);
	j._isMsie = (window.ActiveXObject !== undefined);
  j._isEdge = (navigator.userAgent.indexOf("Edge/") >= 0);
	j._useDataURI = /*!j._isSafari && */ !j._isMsie && !j._isEdge; // safari may be OK here -- untested

  window.requestAnimationFrame || (window.requestAnimationFrame = window.setTimeout);
	for(var i in Jmol) j[i] = Jmol[i]; // allows pre-definition
	return j;
})(document, Jmol);


(function (Jmol, $) {

  Jmol.clazzAlert = function(msg) {
	  //TODO
  }
  Jmol.__$ = $; // local jQuery object -- important if any other module needs to access it (JSmolMenu, for example)

// this library is organized into the following sections:

	// jQuery interface 
	// protected variables
	// feature detection
	// AJAX-related core functionality
	// applet start-up functionality
	// misc core functionality
	// mouse events


	////////////////////// jQuery interface ///////////////////////

	// hooks to jQuery -- if you have a different AJAX tool, feel free to adapt.
	// There should be no other references to jQuery in all the JSmol libraries.

	// automatically switch to returning HTML after the page is loaded
	$(document).ready(function(){ Jmol._document = null });

	Jmol.$ = function(objectOrId, subdiv) {
		// if a subdiv, then return $("#objectOrId._id_subdiv") 
		// or if no subdiv, then just $(objectOrId)
		if (objectOrId == null)alert (subdiv + arguments.callee.caller.toString());
			return $(subdiv ? "#" + objectOrId._id + "_" + subdiv : objectOrId);
	} 

	Jmol._$ = function(id) {
		// either the object or $("#" + id)
		return (typeof id == "string" ? $("#" + id) : id);
	}

	/// special functions:

	Jmol.$ajax = function(info) {
		Jmol._ajaxCall = info.url;
		info.cache = (info.cache != "NO");
    info.url = Jmol._fixProtocol(info.url);
		// don't let jQuery add $_=... to URL unless we 
		// use cache:"NO"; other packages might use $.ajaxSetup() to set this to cache:false
		return $.ajax(info);
	}

  Jmol._fixProtocol = function(url) {
    if (url.indexOf("get3d=True") >= 0)
      url = url.replace(/get3d\=True/, "get3d=true"); // NCI/CADD change 08/2016
  	return (    
    url.indexOf("http://www.rcsb.org/pdb/files/") == 0 && url.indexOf("/ligand/") < 0 ? 
      "http://files.rcsb.org/view/" + url.substring(30).replace(/\.gz/,"")    
    : url.indexOf("http://") == 0 && (
      Jmol._httpProto == "https://"
      || url.indexOf(".gov/") > 0 
      || url.indexOf("http://www.materialsproject") == 0) 
      ? "https" + url.substring(4) : url);
  }
  
	Jmol._getNCIInfo = function(identifier, what, fCallback) {
		return Jmol._getFileData("https://cactus.nci.nih.gov/chemical/structure/"+identifier +"/" + (what == "name" ? "names" : what));
	}
	
	Jmol.$appEvent = function(app, subdiv, evt, f) {
		var o = Jmol.$(app, subdiv); 
		o.off(evt) && f && o.on(evt, f);
	}   

	Jmol.$resize = function (f) {
		return $(window).resize(f);
	}

	//// full identifier expected (could be "body", for example):

	Jmol.$after = function (what, s) {
		return $(what).after(s);
	}

	Jmol.$append = function (what, s) {
		return $(what).append(s);
	}

	Jmol.$bind = function(what, list, f) {
		return (f ? $(what).bind(list, f) : $(what).unbind(list));
	}

	Jmol.$closest = function(what, d) {
		return $(what).closest(d);
	}
	
	Jmol.$get = function(what, i) {
	return $(what).get(i);
	}
 
	// element id expected
			 
	Jmol.$documentOff = function(evt, id) {
		return $(document).off(evt, "#" + id);
	}

	Jmol.$documentOn = function(evt, id, f) {
		return $(document).on(evt, "#" + id, f);
	}

	Jmol.$getAncestorDiv = function(id, className) {
		return $("div." + className + ":has(#" + id + ")")[0];
	}

	Jmol.$supportsIECrossDomainScripting = function() {
		return $.support.iecors;
	}

	//// element id or jQuery object expected:

	Jmol.$attr = function (id, a, val) {
		return Jmol._$(id).attr(a, val);
	}

	Jmol.$css = function(id, style) {
		return Jmol._$(id).css(style);
	}
	 
	Jmol.$find = function(id, d) {
		return Jmol._$(id).find(d);
	}
	
	Jmol.$focus = function(id) {
		return Jmol._$(id).focus();
	}

	Jmol.$html = function(id, html) {
		return Jmol._$(id).html(html);
	}
	 
	Jmol.$offset = function(id) {
		return Jmol._$(id).offset();
	}

	Jmol.$windowOn = function(evt, f) {
		return $(window).on(evt, f);
	}

	Jmol.$prop = function(id, p, val) {
		var o = Jmol._$(id);
		return (arguments.length == 3 ? o.prop(p, val) : o.prop(p));
	}

	Jmol.$remove = function(id) {
		return Jmol._$(id).remove();
	}

	Jmol.$scrollTo = function (id, n) {
		var o = Jmol._$(id);
		return o.scrollTop(n < 0 ? o[0].scrollHeight : n);
	}

	Jmol.$setEnabled = function(id, b) {
		return Jmol._$(id).attr("disabled", b ? null : "disabled");  
	}

  Jmol.$getSize = function(id) {
		var o = Jmol._$(id);
    return [ o.width(), o.height() ]
  }
  
	Jmol.$setSize = function(id, w, h) {
		return Jmol._$(id).width(w).height(h);
	}

  Jmol.$is = function(id, test) { // e.g. :visible
    return Jmol._$(id).is(test);
  }
  
	Jmol.$setVisible = function(id, b) {
		var o = Jmol._$(id);
		return (b ? o.show() : o.hide());  
	}

	Jmol.$submit = function(id) {
		return Jmol._$(id).submit();
	}

	Jmol.$val = function (id, v) {
		var o = Jmol._$(id);
		return (arguments.length == 1 ? o.val() : o.val(v));
	}

	////////////// protected variables ///////////


	Jmol._clearVars = function() {
		// only on page closing -- appears to improve garbage collection

		delete jQuery;
		delete $;
		delete Jmol;
		delete SwingController;
		delete J;
		delete JM;
		delete JS;
		delete JSV;
		delete JU;
		delete JV;
		delete java;
		delete javajs;
		delete Clazz;
		delete c$; // used in p0p; could be gotten rid of
	}

	////////////// feature detection ///////////////

	Jmol.featureDetection = (function(document, window) {

		var features = {};
		features.ua = navigator.userAgent.toLowerCase()

		features.os = (function(){
			var osList = ["linux","unix","mac","win"]
			var i = osList.length;

			while (i--){
				if (features.ua.indexOf(osList[i])!=-1) return osList[i]
			}
			return "unknown";
		})();

		features.browser = function(){
			var ua = features.ua;
			var browserList = ["konqueror","webkit","omniweb","opera","webtv","icab","msie","mozilla"];
			for (var i = 0; i < browserList.length; i++)
			if (ua.indexOf(browserList[i])>=0) 
				return browserList[i];
			return "unknown";
		}
		features.browserName = features.browser();
		features.browserVersion= parseFloat(features.ua.substring(features.ua.indexOf(features.browserName)+features.browserName.length+1));
		features.supportsXhr2 = function() {return ($.support.cors || $.support.iecors)}
		features.allowDestroy = (features.browserName != "msie");
		features.allowHTML5 = (features.browserName != "msie" || navigator.appVersion.indexOf("MSIE 8") < 0);
		features.getDefaultLanguage = function() {
			return navigator.language || navigator.userLanguage || "en-US";
		};

		features._webGLtest = 0;

		features.supportsWebGL = function() {
		if (!Jmol.featureDetection._webGLtest) { 
			var canvas;
			Jmol.featureDetection._webGLtest = ( 
				window.WebGLRenderingContext 
					&& ((canvas = document.createElement("canvas")).getContext("webgl") 
				|| canvas.getContext("experimental-webgl")) ? 1 : -1);
		}
		return (Jmol.featureDetection._webGLtest > 0);
	};

	features.supportsLocalization = function() {
		//<meta charset="utf-8">                                     
		var metas = document.getElementsByTagName('meta'); 
		for (var i= metas.length; --i >= 0;) 
			if (metas[i].outerHTML.toLowerCase().indexOf("utf-8") >= 0) return true;
		return false;
		};

	features.supportsJava = function() {
		if (!Jmol.featureDetection._javaEnabled) {
			if (Jmol._isMsie) {
				if (!navigator.javaEnabled()) {
					Jmol.featureDetection._javaEnabled = -1;
				} else {
					//more likely -- would take huge testing
					Jmol.featureDetection._javaEnabled = 1;
				}
			} else {
				Jmol.featureDetection._javaEnabled = (navigator.javaEnabled() && (!navigator.mimeTypes || navigator.mimeTypes["application/x-java-applet"]) ? 1 : -1);
			}
		}
		return (Jmol.featureDetection._javaEnabled > 0);
	};

	features.compliantBrowser = function() {
		var a = !!document.getElementById;
		var os = features.os;
		// known exceptions (old browsers):
		if (features.browserName == "opera" && features.browserVersion <= 7.54 && os == "mac"
			|| features.browserName == "webkit" && features.browserVersion < 125.12
			|| features.browserName == "msie" && os == "mac"
			|| features.browserName == "konqueror" && features.browserVersion <= 3.3
		) a = false;
		return a;
	}

	features.isFullyCompliant = function() {
		return features.compliantBrowser() && features.supportsJava();
	}

	features.useIEObject = (features.os == "win" && features.browserName == "msie" && features.browserVersion >= 5.5);
	features.useHtml4Object = (features.browserName == "mozilla" && features.browserVersion >= 5) ||
		(features.browserName == "opera" && features.browserVersion >= 8) ||
		(features.browserName == "webkit"/* && features.browserVersion >= 412.2 && features.browserVersion < 500*/); // 500 is a guess; required for 536.3

	features.hasFileReader = (window.File && window.FileReader);

	return features;

})(document, window);


		////////////// AJAX-related core functionality //////////////

	Jmol._ajax = function(info) {
		if (!info.async) {
			return Jmol.$ajax(info).responseText;
		}
		Jmol._ajaxQueue.push(info)
		if (Jmol._ajaxQueue.length == 1)
			Jmol._ajaxDone()
	}
	Jmol._ajaxDone = function() {
		var info = Jmol._ajaxQueue.shift();
		info && Jmol.$ajax(info);
	}

	Jmol._grabberOptions = [
		["$", "NCI(small molecules)"],
		[":", "PubChem(small molecules)"],
		["=", "RCSB(macromolecules)"],
		["*", "PDBe(macromolecules)"]
	];

	Jmol._getGrabberOptions = function(applet) {
		// feel free to adjust this look to anything you want
		if (Jmol._grabberOptions.length == 0)
			return ""


		var s = '<input type="text" id="ID_query" onfocus="jQuery(this).select()" onkeypress="if(13==event.which){Jmol._applets[\'ID\']._search();return false}" size="32" value="" />';
		var b = '<button id="ID_submit" onclick="Jmol._applets[\'ID\']._search()">Search</button></nobr>'
		if (Jmol._grabberOptions.length == 1) {
			s = '<nobr>' + s + '<span style="display:none">';
			b = '</span>' + b;
		} else {
			s += '<br /><nobr>'
		}
		s += '<select id="ID_select">'
		for (var i = 0; i < Jmol._grabberOptions.length; i++) {
			var opt = Jmol._grabberOptions[i];
			s += '<option value="' + opt[0] + '" ' + (i == 0 ? 'selected' : '') + '>' + opt[1] + '</option>';
		}
		s = (s + '</select>' + b).replace(/ID/g, applet._id);
		return '<br />' + s;
	}

	Jmol._getScriptForDatabase = function(database) {
		return (database == "$" ? Jmol.db._nciLoadScript : database == ":" ? Jmol.db._pubChemLoadScript : Jmol.db._fileLoadScript);
	}

	 //   <dataset><record><structureId>1BLU</structureId><structureTitle>STRUCTURE OF THE 2[4FE-4S] FERREDOXIN FROM CHROMATIUM VINOSUM</structureTitle></record><record><structureId>3EUN</structureId><structureTitle>Crystal structure of the 2[4Fe-4S] C57A ferredoxin variant from allochromatium vinosum</structureTitle></record></dataset>

	Jmol._setInfo = function(applet, database, data) {
		var info = [];
		var header = "";
		if (data.indexOf("ERROR") == 0)
			header = data;
		else
			switch (database) {
			case "=":
				var S = data.split("<dimStructure.structureId>");
				var info = ["<table>"];
				for (var i = 1; i < S.length; i++) {
					info.push("<tr><td valign=top><a href=\"javascript:Jmol.search(" + applet._id + ",'=" + S[i].substring(0, 4) + "')\">" + S[i].substring(0, 4) + "</a></td>");
					info.push("<td>" + S[i].split("Title>")[1].split("</")[0] + "</td></tr>");
				}
				info.push("</table>");
				header = (S.length - 1) + " matches";
				break;      
			case "$": // NCI
			case ":": // pubChem
			break;
			default:
				return;
		}
		applet._infoHeader = header;
		applet._info = info.join("");
		applet._showInfo(true);
	}

	Jmol._loadSuccess = function(a, fSuccess) {
		if (!fSuccess)
			return;
		Jmol._ajaxDone();
		fSuccess(a);
	}

	Jmol._loadError = function(fError){
		Jmol._ajaxDone();
		Jmol.say("Error connecting to server: " + Jmol._ajaxCall);  
		null!=fError&&fError()
	}

	Jmol._isDatabaseCall = function(query) {
		return (Jmol.db._databasePrefixes.indexOf(query.substring(0, 1)) >= 0);
	}

	Jmol._getDirectDatabaseCall = function(query, checkXhr2) {
		if (checkXhr2 && !Jmol.featureDetection.supportsXhr2())
			return query;
		var pt = 2;
		var db = query.substring(0,pt)
		var call = Jmol.db._DirectDatabaseCalls[db] || Jmol.db._DirectDatabaseCalls[db = query.substring(0,--pt)];
		if (call) {
      // one of the special set :, =, $, ==
			if (db == ":") {
        // PubChem
				var ql = query.toLowerCase();
				if (!isNaN(parseInt(query.substring(1)))) {
					query = "cid/" + query.substring(1);
				} else if (ql.indexOf(":smiles:") == 0) {
					call += "?POST?smiles=" + query.substring(8);
					query = "smiles";
				} else if (ql.indexOf(":cid:") == 0) {
					query = "cid/" + query.substring(5);
				} else {
					if (ql.indexOf(":name:") == 0)
						query = query.substring(5);
					else if (ql.indexOf(":cas:") == 0)
						query = query.substring(4);
					query = "name/" + encodeURIComponent(query.substring(pt));
				}
			} else {
				query = encodeURIComponent(query.substring(pt));		
			}
      if (query.indexOf(".mmtf") >= 0) {
        query = "https://mmtf.rcsb.org/v1.0/full/" + query.replace(/\.mmtf/, "");
			} else if (call.indexOf("FILENCI") >= 0) {
				query = query.replace(/\%2F/g, "/");				
				query = call.replace(/\%FILENCI/, query);
			} else {
				query = call.replace(/\%FILE/, query);
			}
		}
		return query;
	}

	Jmol._getRawDataFromServer = function(database,query,fSuccess,fError,asBase64,noScript){
	  // note that this method is now only enabled for "_"
	  // server-side processing of database queries was too slow and only useful for 
	  // the IMAGE option, which has been abandoned.
		var s = 
			"?call=getRawDataFromDatabase&database=" + database + (query.indexOf("?POST?") >= 0 ? "?POST?" : "")
				+ "&query=" + encodeURIComponent(query)
				+ (asBase64 ? "&encoding=base64" : "")
				+ (noScript ? "" : "&script=" + encodeURIComponent(Jmol._getScriptForDatabase(database)));
		return Jmol._contactServer(s, fSuccess, fError);
	}

	Jmol._checkFileName = function(applet, fileName, isRawRet) {
		if (Jmol._isDatabaseCall(fileName)) {
			if (isRawRet)
				Jmol._setQueryTerm(applet, fileName);
			fileName = Jmol._getDirectDatabaseCall(fileName, true);
			if (Jmol._isDatabaseCall(fileName)) {
				// xhr2 not supported (MSIE)
				fileName = Jmol._getDirectDatabaseCall(fileName, false);
				isRawRet && (isRawRet[0] = true);
			}
		}
		return fileName;
	}
	
	Jmol._checkCache = function(applet, fileName, fSuccess) {
		if (applet._cacheFiles && Jmol._fileCache && !fileName.endsWith(".js")) {
			var data = Jmol._fileCache[fileName];
			if (data) {
				System.out.println("using "  + data.length + " bytes of cached data for "  + fileName);
				fSuccess(data);
				return null;
			} else {
				fSuccess = function(fileName, data) { fSuccess(Jmol._fileCache[fileName] = data) };     
			}
		}
		return fSuccess;
	}
	
  Jmol.playAudio = function(filePath) {
    Jmol.playAudio(null, filePath);
  }
  
  Jmol.playAudio = function(applet, params) {
  
    var get = (params.get ? function(key){return params.get(key)} : null);
    var put = (params.put ? function(key,val){return params.put(key,val)} : null);
    var filePath = (get ? get("audioFile") : params);
    var jmolAudio = get && get("audioPlayer");
    var audio = document.createElement("audio");
    put && put("audioElement", audio);
    var callback = null;
    if (jmolAudio) {
      callback = function(type){jmolAudio.processUpdate(type)};
      jmolAudio.myClip = {
         open: function() {callback("open")},
         start: function() { audio.play(); callback("start")},
         loop: function(n) { audio.loop = (n != 0); },
         stop: function() { audio.pause(); },
         close: function() { callback("close") },
         setMicrosecondPosition: function(us) { audio.currentTime = us / 1e6; }
      };
    }    
    audio.controls = "true";
    audio.src = filePath;
    if (get && get("loop"))
      audio.loop = "true";
    if (callback) {
      audio.addEventListener("pause", function() {
          callback("pause");
      });
      audio.addEventListener("play", function() {
          callback("play");
      });
      audio.addEventListener("playing", function() {
          callback("playing");
      });
      audio.addEventListener("ended", function() {
          callback("ended");
      });
      callback("open")
    }
  }
  
	Jmol._loadFileData = function(applet, fileName, fSuccess, fError){
		var isRaw = [];
		fileName = Jmol._checkFileName(applet, fileName, isRaw);
		fSuccess = Jmol._checkCache(applet, fileName, fSuccess);
		if (isRaw[0]) {
				Jmol._getRawDataFromServer("_",fileName,fSuccess,fError);   
				return;
		} 
		var info = {
			type: "GET",
			dataType: "text",
			url: fileName,
			async: Jmol._asynchronous,
			success: function(a) {Jmol._loadSuccess(a, fSuccess)},
			error: function() {Jmol._loadError(fError)}
		}
		Jmol._checkAjaxPost(info);
		Jmol._ajax(info);
	}

	Jmol._getInfoFromDatabase = function(applet, database, query){
		if (database == "====") {
			var data = Jmol.db._restQueryXml.replace(/QUERY/,query);
			var info = {
				dataType: "text",
				type: "POST",
				contentType:"application/x-www-form-urlencoded",
				url: Jmol.db._restQueryUrl,
				data: encodeURIComponent(data) + "&req=browser",
				success: function(data) {Jmol._ajaxDone();Jmol._extractInfoFromRCSB(applet, database, query, data)},
				error: function() {Jmol._loadError(null)},
				async: Jmol._asynchronous
			}
			return Jmol._ajax(info);
		}   
		query = "?call=getInfoFromDatabase&database=" + database
				+ "&query=" + encodeURIComponent(query);
		return Jmol._contactServer(query, function(data) {Jmol._setInfo(applet, database, data)});
	}

	Jmol._extractInfoFromRCSB = function(applet, database, query, output) {
		var n = output.length/5;
		if (n == 0)
			return; 
		if (query.length == 4 && n != 1) {
			var QQQQ = query.toUpperCase();
			var pt = output.indexOf(QQQQ);
			if (pt > 0 && "123456789".indexOf(QQQQ.substring(0, 1)) >= 0)
				output = QQQQ + "," + output.substring(0, pt) + output.substring(pt + 5);
			if (n > 50)
				output = output.substring(0, 250);
			output = output.replace(/\n/g,",");
			var url = Jmol._restReportUrl.replace(/IDLIST/,output);
			Jmol._loadFileData(applet, url, function(data) {Jmol._setInfo(applet, database, data) });   
		}
	}

	Jmol._checkAjaxPost = function(info) {
		var pt = info.url.indexOf("?POST?");
		if (pt > 0) {
			info.data = info.url.substring(pt + 6);
			info.url = info.url.substring(0, pt);
			info.type = "POST";
			info.contentType = "application/x-www-form-urlencoded";
		}
	}
	Jmol._contactServer = function(data,fSuccess,fError){
		var info = {
			dataType: "text",
			type: "GET",
			url: Jmol._serverUrl + data,
			success: function(a) {Jmol._loadSuccess(a, fSuccess)},
			error:function() { Jmol._loadError(fError) },
			async:fSuccess ? Jmol._asynchronous : false
		}
		Jmol._checkAjaxPost(info);
		return Jmol._ajax(info);
	}

	Jmol._setQueryTerm = function(applet, query) {
		if (!query || !applet._hasOptions || query.substring(0, 7) == "http://")
			return;
		if (Jmol._isDatabaseCall(query)) {
			var database = query.substring(0, 1);
			query = query.substring(1);
			if (query.substring(0,1) == database && "=$".indexOf(database) >= 0)
				query = query.substring(1);
			var d = Jmol._getElement(applet, "select");
			if (d && d.options)
				for (var i = 0; i < d.options.length; i++)
					if (d[i].value == database)
						d[i].selected = true;
		}
		Jmol.$val(Jmol.$(applet, "query"), query);
	}

	Jmol._search = function(applet, query, script) {
		arguments.length > 1 || (query = null);
		Jmol._setQueryTerm(applet, query);
		query || (query = Jmol.$val(Jmol.$(applet, "query")))
		if (query.indexOf("!") == 0) {
		// command prompt in this box as well
		// remove exclamation point "immediate" indicator
			applet._script(query.substring(1));
			return;
		} 
		query && (query = query.replace(/\"/g, ""));
		applet._showInfo(false);
		Jmol._searchMol(applet, query, script, true);
	}

	Jmol._searchMol = function(applet, query, script, checkView) {
		var database;
		if (Jmol._isDatabaseCall(query)) {
			database = query.substring(0, 1);
			query = query.substring(1);
		} else {
			database = (applet._hasOptions ? Jmol.$val(Jmol.$(applet, "select")) : "$");
		}
		if (database == "=" && query.length == 3)
			query = "=" + query; // this is a ligand      
		var dm = database + query;
		if (!query || dm.indexOf("?") < 0 && dm == applet._thisJmolModel) {
			return;    
		}
		applet._thisJmolModel = dm;
		var view;
		if (checkView && applet._viewSet != null && (view = Jmol.View.__findView(applet._viewSet, {chemID:dm})) != null) {
			Jmol.View.__setView(view, applet, false);
			return;
		}

		if (database == "$" || database == ":")
			applet._jmolFileType = "MOL";
		else if (database == "=")
			applet._jmolFileType = "PDB";
		applet._searchDatabase(query, database, script);
	}

	Jmol._searchDatabase = function(applet, query, database, script) {
		applet._showInfo(false);
		if (query.indexOf("?") >= 0) {
			Jmol._getInfoFromDatabase(applet, database, query.split("?")[0]);
			return true;
		}
		if (Jmol.db._DirectDatabaseCalls[database]) {
			applet._loadFile(database + query, script);
			return true;
		}
		return false;
	}

	Jmol._syncBinaryOK="?";

	Jmol._canSyncBinary = function(isSilent) {
		if (Jmol._isAsync) return true;
		if (self.VBArray) return (Jmol._syncBinaryOK = false);
		if (Jmol._syncBinaryOK != "?") return Jmol._syncBinaryOK;
		Jmol._syncBinaryOK = true;
		try {
			var xhr = new window.XMLHttpRequest();
			xhr.open( "text", Jmol._ajaxTestSite, false );
			if (xhr.hasOwnProperty("responseType")) {
				xhr.responseType = "arraybuffer";
			} else if (xhr.overrideMimeType) {
				xhr.overrideMimeType('text/plain; charset=x-user-defined');
			}
		} catch( e ) {
			var s = "JSmolCore.js: synchronous binary file transfer is requested but not available";
			System.out.println(s);
			if (Jmol._alertNoBinary && !isSilent)
				alert (s)
			return Jmol._syncBinaryOK = false;
		}
		return true;  
	}

	Jmol._binaryTypes = ["mmtf",".gz",".bz2",".jpg",".gif",".png",".zip",".jmol",".bin",".smol",".spartan",".pmb",".mrc",".map",".ccp4",".dn6",".delphi",".omap",".pse",".dcd",".lut",".uk/pdbe/densities/"];

	Jmol.isBinaryUrl = function(url) {
		for (var i = Jmol._binaryTypes.length; --i >= 0;)
			if (url.indexOf(Jmol._binaryTypes[i]) >= 0) return true;
		return false;
	}

	Jmol._getFileData = function(fileName, fSuccess, doProcess) {
		// use host-server PHP relay if not from this host
		var isBinary = Jmol.isBinaryUrl(fileName);
		var isPDB = (fileName.indexOf(".gz") >= 0 && fileName.indexOf("rcsb.org") >= 0);
		if (isPDB) {
			// avoid unnecessary binary transfer
			fileName = fileName.replace(/\.gz/,"");
			isBinary = false;
		}
		var asBase64 = (isBinary && !fSuccess && !Jmol._canSyncBinary(isPDB));
		var isPost = (fileName.indexOf("?POST?") >= 0);
		if (fileName.indexOf("file:/") == 0 && fileName.indexOf("file:///") != 0)
			fileName = "file://" + fileName.substring(5);      /// fixes IE problem
		var isMyHost = (fileName.indexOf("://") < 0 || fileName.indexOf(document.location.protocol) == 0 && fileName.indexOf(document.location.host) >= 0);
    var isHttps2Http = (Jmol._httpProto == "https://" && fileName.indexOf("http://") == 0);
		var isDirectCall = Jmol._isDirectCall(fileName);
    if (!isDirectCall && fileName.indexOf("?ALLOWSORIGIN?") >= 0) {
      isDirectCall = true;
			fileName = fileName.replace(/\?ALLOWSORIGIN\?/,"");
    }
		//if (fileName.indexOf("http://pubchem.ncbi.nlm.nih.gov/") == 0)isDirectCall = false;

		var cantDoSynchronousLoad = (!isMyHost && Jmol.$supportsIECrossDomainScripting());
		var data = null;
		if (isHttps2Http || asBase64 || !isMyHost && !isDirectCall || !fSuccess && cantDoSynchronousLoad ) {
				data = Jmol._getRawDataFromServer("_",fileName, fSuccess, fSuccess, asBase64, true);
		} else {
			fileName = fileName.replace(/file:\/\/\/\//, "file://"); // opera
			var info = {dataType:(isBinary ? "binary" : "text"),async:!!fSuccess};
			if (isPost) {
				info.type = "POST";
				info.url = fileName.split("?POST?")[0]
				info.data = fileName.split("?POST?")[1]
			} else {
				info.type = "GET";
				info.url = fileName;
			}
			if (fSuccess) {
				info.success = function(data) { fSuccess(Jmol._xhrReturn(info.xhr))};
				info.error = function() { fSuccess(info.xhr.statusText)};
			}
			info.xhr = Jmol.$ajax(info);
			if (!fSuccess) {
				data = Jmol._xhrReturn(info.xhr);
			}
		}
		if (!doProcess)
			return data;
		if (data == null) {
			data = "";
			isBinary = false;
		}
		isBinary && (isBinary = Jmol._canSyncBinary(true));
		return (isBinary ? Jmol._strToBytes(data) : JU.SB.newS(data));
	}
	
	Jmol._xhrReturn = function(xhr){
		if (!xhr.responseText && xhr.responseText !== '' || self.Clazz && Clazz.instanceOf(xhr.response, self.ArrayBuffer)) {
			// Safari or error 
			return xhr.response || xhr.statusText;
		} 
		return xhr.responseText;
	}

	Jmol._isDirectCall = function(url) {
    if (url.indexOf("?ALLOWSORIGIN?") >= 0)
      return true;
		for (var key in Jmol.db._DirectDatabaseCalls) {
			if (key.indexOf(".") >= 0 && url.indexOf(key) >= 0)
				return true;
		}
		return false;
	}

	Jmol._cleanFileData = function(data) {
		if (data.indexOf("\r") >= 0 && data.indexOf("\n") >= 0) {
			return data.replace(/\r\n/g,"\n");
		}
		if (data.indexOf("\r") >= 0) {
			return data.replace(/\r/g,"\n");
		}
		return data;
	};

	Jmol._getFileType = function(name) {
		var database = name.substring(0, 1);
		if (database == "$" || database == ":")
			return "MOL";
		if (database == "=")
			return (name.substring(1,2) == "=" ? "LCIF" : "PDB");
		// just the extension, which must be PDB, XYZ..., CIF, or MOL
		name = name.split('.').pop().toUpperCase();
		return name.substring(0, Math.min(name.length, 3));
	};

	Jmol._getZ = function(applet, what) {
		return applet && applet._z && applet._z[what] || Jmol._z[what];
	}
	
	Jmol._incrZ = function(applet, what) {
		return applet && applet._z && ++applet._z[what] || ++Jmol._z[what];
	}
	
  Jmol._hideLocalFileReader = function(applet, cursor) {
    applet._localReader && Jmol.$setVisible(applet._localReader, false);
    applet._readingLocal = false;
    Jmol.setCursor(applet, 0);
  }
  
  Jmol.loadFileFromDialog = function(applet) {
    Jmol.loadFileAsynchronously(null, applet, null, null);
  }
  
	Jmol.loadFileAsynchronously = function(fileLoadThread, applet, fileName, appData) {
		if (fileName && fileName.indexOf("?") != 0) {
			// LOAD ASYNC command
			var fileName0 = fileName;
			fileName = Jmol._checkFileName(applet, fileName);
			var fSuccess = function(data) {Jmol._setData(fileLoadThread, fileName, fileName0, data, appData, applet)};
			fSuccess = Jmol._checkCache(applet, fileName, fSuccess);
			if (fileName.indexOf("|") >= 0)
				fileName = fileName.split("|")[0];
			return (fSuccess == null ? null : Jmol._getFileData(fileName, fSuccess));		
		}
		// we actually cannot suggest a fileName, I believe.
		if (!Jmol.featureDetection.hasFileReader) {
        var msg = "Local file reading is not enabled in your browser";
				return (fileLoadThread ? fileLoadThread.setData(msg, null, null, appData, applet) : alert(msg));
    }
		if (!applet._localReader) {
			var div = '<div id="ID" style="z-index:'+Jmol._getZ(applet, "fileOpener") + ';position:absolute;background:#E0E0E0;left:10px;top:10px"><div style="margin:5px 5px 5px 5px;"><button id="ID_loadurl">URL</button><input type="file" id="ID_files" /><button id="ID_loadfile">load</button><button id="ID_cancel">cancel</button></div><div>'
			Jmol.$after("#" + applet._id + "_appletdiv", div.replace(/ID/g, applet._id + "_localReader"));
			applet._localReader = Jmol.$(applet, "localReader");
		}
		Jmol.$appEvent(applet, "localReader_loadurl", "click");
		Jmol.$appEvent(applet, "localReader_loadurl", "click", function(evt) {
      var file = prompt("Enter a URL");
      if (!file)return
      Jmol._hideLocalFileReader(applet, 0);
      Jmol._setData(null, file, file, null, appData, applet);
		});
		Jmol.$appEvent(applet, "localReader_loadfile", "click");
		Jmol.$appEvent(applet, "localReader_loadfile", "click", function(evt) {
			var file = Jmol.$(applet, "localReader_files")[0].files[0];   
			var reader = new FileReader();
			reader.onloadend = function(evt) {
				if (evt.target.readyState == FileReader.DONE) { // DONE == 2
          Jmol._hideLocalFileReader(applet, 0);
					Jmol._setData(fileLoadThread, file.name, file.name, evt.target.result, appData, applet);
				}
			};
      try {
			  reader.readAsArrayBuffer(file);
      } catch (e) {
        alert("You must select a file first.");
      }
		});
		Jmol.$appEvent(applet, "localReader_cancel", "click");
		Jmol.$appEvent(applet, "localReader_cancel", "click", function(evt) {
      Jmol._hideLocalFileReader(applet);
      if (fileLoadThread)
  			fileLoadThread.setData("#CANCELED#", null, null, appData, applet);
		});
		Jmol.$setVisible(applet._localReader, true);
    applet._readingLocal = true;
	}

  Jmol._setData = function(fileLoadThread, filename, filename0, data, appData, applet) {
  	data && (data = Jmol._strToBytes(data));
		if (data != null && (fileLoadThread == null || filename.indexOf(".jdx") >= 0))
			Jmol.Cache.put("cache://" + filename, data); 
    if (fileLoadThread == null) {
      applet._applet.openFileAsyncSpecial(data == null ? filename : "cache://" + filename, 1);
    } else {          
		  fileLoadThread.setData(filename, filename0, data, appData);
    }    
  }

	Jmol.doAjax = function(url, postOut, dataOut) {
		// called by org.jmol.awtjs2d.JmolURLConnection.doAjax()
		url = url.toString();
		if (dataOut) {

			if (url.indexOf("http://") != 0 && url.indexOf("https://") != 0)
				return Jmol._saveFile(url, dataOut);
			var info = {async:false,url:url,type:"POST",
				data:typeof data == "string" ? dataOut :";base64," + (JU || J.util).Base64.getBase64(dataOut).toString(), processData:false
			};
			return Jmol.$ajax(info).responseText;
		}
		if (postOut)
			url += "?POST?" + postOut;
		return Jmol._getFileData(url, null, true);
	}

	// Jmol._localFileSaveFunction --  // do something local here; Maybe try the FileSave interface? return true if successful
	 
	Jmol._saveFile = function(filename, data, mimetype, encoding) {
		if (Jmol._localFileSaveFunction && Jmol._localFileSaveFunction(filename, data))
			return "OK";
		var filename = filename.substring(filename.lastIndexOf("/") + 1);
		mimetype || (mimetype = (filename.indexOf(".pdf") >= 0 ? "application/pdf" 
			: filename.indexOf(".png") >= 0 ? "image/png" 
			: filename.indexOf(".gif") >= 0 ? "image/gif" 
			: filename.indexOf(".jpg") >= 0 ? "image/jpg" : ""));
		var isString = (typeof data == "string");
		
	if (isString && data.indexOf(";base64,") >= 0) {
	  data = data.split(";base64,")[1];
	} else {	
   	data = (JU || J.util).Base64.getBase64(isString ? data.getBytes("UTF-8") : data).toString();
   	}
		encoding || (encoding = "base64");
		
		var url = Jmol._serverUrl;
		url && url.indexOf("your.server") >= 0 && (url = "");
		if (Jmol._useDataURI || !url) {
			// Asynchronous output generated using an anchor tag
			// btoa does not work with UTF-8 data///encoding || (data = btoa(data));
			var a = document.createElement("a");
			a.href = "data:" + mimetype + ";base64," + data;
			a.type = mimetype || (mimetype = "text/plain;charset=utf-8");	
			a.download = filename;
			a.target = "_blank";
				$("body").append(a);
			a.click();
			a.remove();		
		} else {
		// Asynchronous outputto be reflected as a download
			if (!Jmol._formdiv) {
					var sform = '<div id="__jsmolformdiv__" style="display:none">\
						<form id="__jsmolform__" method="post" target="_blank" action="">\
						<input name="call" value="saveFile"/>\
						<input id="__jsmolmimetype__" name="mimetype" value=""/>\
						<input id="__jsmolencoding__" name="encoding" value=""/>\
						<input id="__jsmolfilename__" name="filename" value=""/>\
						<textarea id="__jsmoldata__" name="data"></textarea>\
						</form>\
						</div>'
				Jmol.$after("body", sform);
				Jmol._formdiv = "__jsmolform__";
			}
			Jmol.$attr(Jmol._formdiv, "action", url + "?" + (new Date()).getMilliseconds());
			Jmol.$val("__jsmoldata__", data);
			Jmol.$val("__jsmolfilename__", filename);
			Jmol.$val("__jsmolmimetype__", mimetype);
			Jmol.$val("__jsmolencoding__", encoding);
			Jmol.$submit("__jsmolform__");
			Jmol.$val("__jsmoldata__", "");
			Jmol.$val("__jsmolfilename__", "");
		}
		return "OK";
	}

	Jmol._strToBytes = function(s) {
		if (Clazz.instanceOf(s, self.ArrayBuffer))
			return Clazz.newByteArray(-1, s);
		var b = Clazz.newByteArray(s.length, 0);
		for (var i = s.length; --i >= 0;)
			b[i] = s.charCodeAt(i) & 0xFF;
		return b;
	}

	////////////// applet start-up functionality //////////////

	Jmol._setConsoleDiv = function (d) {
		if (!self.Clazz)return;
		Clazz.setConsoleDiv(d);
	}

	Jmol._registerApplet = function(id, applet) {
		return window[id] = Jmol._applets[id] = Jmol._applets.master = Jmol._applets[id + "__" + Jmol._syncId + "__"] = applet;
	} 

	Jmol._readyCallback = function (appId,fullId,isReady,javaApplet,javaAppletPanel) {
		appId = appId.split("_object")[0];
    var applet = Jmol._applets[appId];
		isReady = (isReady.booleanValue ? isReady.booleanValue() : isReady);
		// necessary for MSIE in strict mode -- apparently, we can't call 
		// jmol._readyCallback, but we can call Jmol._readyCallback. Go figure...
    if (isReady) {
      // when leaving page, Java applet may be dead 
      applet._appletPanel = (javaAppletPanel || javaApplet);
      applet._applet = javaApplet;
    }
		Jmol._track(applet)._readyCallback(appId, fullId, isReady);
	}

	Jmol._getWrapper = function(applet, isHeader) {

			// id_appletinfotablediv
			//     id_appletdiv
			//     id_coverdiv
			//     id_infotablediv
			//       id_infoheaderdiv
			//          id_infoheaderspan
			//          id_infocheckboxspan
			//       id_infodiv


			// for whatever reason, without DOCTYPE, with MSIE, "height:auto" does not work, 
			// and the text scrolls off the page.
			// So I'm using height:95% here.
			// The table was a fix for MSIE with no DOCTYPE tag to fix the miscalculation
			// in height of the div when using 95% for height. 
			// But it turns out the table has problems with DOCTYPE tags, so that's out. 
			// The 95% is a compromise that we need until the no-DOCTYPE MSIE solution is found. 
			// (100% does not work with the JME linked applet)
		var s;
		// ... here are just for clarification in this code; they are removed immediately
		if (isHeader) {
			var img = ""; 
			if (applet._coverImage){
				var more = " onclick=\"Jmol.coverApplet(ID, false)\" title=\"" + (applet._coverTitle) + "\"";
				var play = "<image id=\"ID_coverclickgo\" src=\"" + applet._makeLiveImage + "\" style=\"width:25px;height:25px;position:absolute;bottom:10px;left:10px;"
					+ "z-index:" + Jmol._getZ(applet, "coverImage")+";opacity:0.5;\"" + more + " />"  
				img = "<div id=\"ID_coverdiv\" style=\"background-color:red;z-index:" + Jmol._getZ(applet, "coverImage")+";width:100%;height:100%;display:inline;position:absolute;top:0px;left:0px\"><image id=\"ID_coverimage\" src=\""
				 + applet._coverImage + "\" style=\"width:100%;height:100%\"" + more + "/>" + play + "</div>";
			}
      
      var wait = (applet._isJava ? "" : "<image id=\"ID_waitimage\" src=\"" + applet._j2sPath + "/img/cursor_wait.gif\" style=\"display:none;position:absolute;bottom:10px;left:10px;"
					+ "z-index:" + Jmol._getZ(applet, "coverImage")+";\" />");  

			var css = Jmol._appletCssText.replace(/\'/g,'"');
      var spinner = applet._getSpinner && applet._getSpinner();
      applet._spinner = spinner = (!spinner || spinner == "none" ? "" : "background-image:url("+spinner + "); background-repeat:no-repeat; background-position:center;");    
			css = spinner + (css.indexOf("style=\"") >= 0 ? css.split("style=\"")[1] : "\" " + css);
			s = "\
...<div id=\"ID_appletinfotablediv\" style=\"width:Wpx;height:Hpx;position:relative;font-size:14px;text-align:left\">IMG WAIT\
......<div id=\"ID_appletdiv\" style=\"z-index:" + Jmol._getZ(applet, "header") + ";width:100%;height:100%;position:absolute;top:0px;left:0px;" + css + ">";
			var height = applet._height;
			var width = applet._width;
			if (typeof height !== "string" || height.indexOf("%") < 0) 
				height += "px";
			if (typeof width !== "string" || width.indexOf("%") < 0)
				width += "px";
			s = s.replace(/IMG/, img).replace(/WAIT/, wait).replace(/Hpx/g, height).replace(/Wpx/g, width);
		} else {
			s = "\
......</div>\
......<div id=\"ID_2dappletdiv\" style=\"position:absolute;width:100%;height:100%;overflow:hidden;display:none\"></div>\
......<div id=\"ID_infotablediv\" style=\"width:100%;height:100%;position:absolute;top:0px;left:0px\">\
.........<div id=\"ID_infoheaderdiv\" style=\"height:20px;width:100%;background:yellow;display:none\"><span id=\"ID_infoheaderspan\"></span><span id=\"ID_infocheckboxspan\" style=\"position:absolute;text-align:right;right:1px;\"><a href=\"javascript:Jmol.showInfo(ID,false)\">[x]</a></span></div>\
.........<div id=\"ID_infodiv\" style=\"position:absolute;top:20px;bottom:0px;width:100%;height:100%;overflow:auto\"></div>\
......</div>\
...</div>";
		}
		return s.replace(/\.\.\./g,"").replace(/[\n\r]/g,"").replace(/ID/g, applet._id);
	}

  Jmol._hideLoadingSpinner = function(applet) {
    if (applet._spinner)
      Jmol.$css(Jmol.$(applet, "appletdiv"), {"background-image": ""});
  }
  
	Jmol._documentWrite = function(text) {
		if (Jmol._document) {
			if (Jmol._isXHTML && !Jmol._XhtmlElement) {
				var s = document.getElementsByTagName("script");
				Jmol._XhtmlElement = s.item(s.length - 1);
				Jmol._XhtmlAppendChild = false;
			}
			if (Jmol._XhtmlElement)
				Jmol._domWrite(text);
			else
				Jmol._document.write(text);
		}
		return text;
	}

	Jmol._domWrite = function(data) {
		var pt = 0
		var Ptr = []
		Ptr[0] = 0
		while (Ptr[0] < data.length) {
			var child = Jmol._getDomElement(data, Ptr);
			if (!child)
				break;
			if (Jmol._XhtmlAppendChild)
				Jmol._XhtmlElement.appendChild(child);
			else
				Jmol._XhtmlElement.parentNode.insertBefore(child, _jmol.XhtmlElement);
		}
	}

	Jmol._getDomElement = function(data, Ptr, closetag, lvel) {

		// there is no "document.write" in XHTML

		var e = document.createElement("span");
		e.innerHTML = data;
		Ptr[0] = data.length;

/*
	// unnecessary ?  

		closetag || (closetag = "");
		lvel || (lvel = 0);
		var pt0 = Ptr[0];
		var pt = pt0;
		while (pt < data.length && data.charAt(pt) != "<") 
			pt++
		if (pt != pt0) {
			var text = data.substring(pt0, pt);
			Ptr[0] = pt;
			return document.createTextNode(text);
		}
		pt0 = ++pt;
		var ch;
		while (pt < data.length && "\n\r\t >".indexOf(ch = data.charAt(pt)) < 0) 
			pt++;
		var tagname = data.substring(pt0, pt);
		var e = (tagname == closetag  || tagname == "/" ? ""
			: document.createElementNS ? document.createElementNS('http://www.w3.org/1999/xhtml', tagname)
			: document.createElement(tagname));
		if (ch == ">") {
			Ptr[0] = ++pt;
			return e;
		}
		while (pt < data.length && (ch = data.charAt(pt)) != ">") {
			while (pt < data.length && "\n\r\t ".indexOf(ch = data.charAt(pt)) >= 0) 
				pt++;
			pt0 = pt;
			while (pt < data.length && "\n\r\t =/>".indexOf(ch = data.charAt(pt)) < 0) 
				pt++;
			var attrname = data.substring(pt0, pt).toLowerCase();
			if (attrname && ch != "=")
				e.setAttribute(attrname, "true");
			while (pt < data.length && "\n\r\t ".indexOf(ch = data.charAt(pt)) >= 0) 
				pt++;
			if (ch == "/") {
				Ptr[0] = pt + 2;
				return e;
			} else if (ch == "=") {
				var quote = data.charAt(++pt);
				pt0 = ++pt;
				while (pt < data.length && (ch = data.charAt(pt)) != quote) 
					pt++;
				var attrvalue = data.substring(pt0, pt);
				e.setAttribute(attrname, attrvalue);
				pt++;
			}
		}
		Ptr[0] = ++pt;
		while (Ptr[0] < data.length) {
			var child = Jmol._getDomElement(data, Ptr, "/" + tagname, lvel+1);
			if (!child)
				break;
			e.appendChild(child);
		}
*/
		return e;    
	}

	Jmol._setObject = function(obj, id, Info) {
		obj._id = id;
		obj.__Info = {};
		Info.z && Info.zIndexBase && (Jmol._z = Jmol._getZOrders(Info.zIndexBase));
		for (var i in Info)
			obj.__Info[i] = Info[i];
		(obj._z = Info.z) || Info.zIndexBase && (obj._z = obj.__Info.z = Jmol._getZOrders(Info.zIndexBase));
		obj._width = Info.width;
		obj._height = Info.height;
		obj._noscript = !obj._isJava && Info.noscript;
		obj._console = Info.console;
		obj._cacheFiles = !!Info.cacheFiles;
		obj._viewSet = (Info.viewSet == null || obj._isJava ? null : "Set" + Info.viewSet);
		if (obj._viewSet != null) {
			Jmol.View.__init(obj);
			obj._currentView = null;
		}
		!Jmol._fileCache && obj._cacheFiles && (Jmol._fileCache = {});
		if (!obj._console)
			obj._console = obj._id + "_infodiv";
		if (obj._console == "none")
			obj._console = null;

		obj._color = (Info.color ? Info.color.replace(/0x/,"#") : "#FFFFFF");
		obj._disableInitialConsole = Info.disableInitialConsole;
		obj._noMonitor = Info.disableJ2SLoadMonitor;
		Jmol._j2sPath && (Info.j2sPath = Jmol._j2sPath);
		obj._j2sPath = Info.j2sPath;
		obj._coverImage = Info.coverImage;
    obj._makeLiveImage = Info.makeLiveImage || Info.j2sPath + "/img/play_make_live.jpg";
		obj._isCovered = !!obj._coverImage; 
		obj._deferApplet = Info.deferApplet || obj._isCovered && obj._isJava; // must do this if covered in Java
		obj._deferUncover = Info.deferUncover && !obj._isJava; // can't do this with Java
		obj._coverScript = Info.coverScript;
		obj._coverTitle = Info.coverTitle;

		if (!obj._coverTitle)
			obj._coverTitle = (obj._deferApplet ? "activate 3D model" : "3D model is loading...")
		obj._containerWidth = obj._width + ((obj._width==parseFloat(obj._width))? "px":"");
		obj._containerHeight = obj._height + ((obj._height==parseFloat(obj._height))? "px":"");
		obj._info = "";
		obj._infoHeader = obj._jmolType + ' "' + obj._id + '"'
		obj._hasOptions = Info.addSelectionOptions;
		obj._defaultModel = Info.defaultModel;
		obj._readyScript = (Info.script ? Info.script : "");
		obj._readyFunction = Info.readyFunction;
		if (obj._coverImage && !obj._deferApplet)
			obj._readyScript += ";javascript " + id + "._displayCoverImage(false)";
		obj._src = Info.src;

	}

	Jmol._addDefaultInfo = function(Info, DefaultInfo) {
		for (var x in DefaultInfo)
			if (typeof Info[x] == "undefined")
				Info[x] = DefaultInfo[x];
		Jmol._use && (Info.use = Jmol._use);
		if (Info.use.indexOf("SIGNED") >= 0) {
			if (Info.jarFile.indexOf("Signed") < 0)
				Info.jarFile = Info.jarFile.replace(/Applet/,"AppletSigned");
			Info.use = Info.use.replace(/SIGNED/, "JAVA");
			Info.isSigned = true;
		}
	}

	Jmol._syncedApplets = [];
	Jmol._syncedCommands = [];
	Jmol._syncedReady = [];
	Jmol._syncReady = false;
	Jmol._isJmolJSVSync = false;

	Jmol._setReady = function(applet) {
		Jmol._syncedReady[applet] = 1;
		var n = 0;
		for (var i = 0; i < Jmol._syncedApplets.length; i++) {
			if (Jmol._syncedApplets[i] == applet._id) {
				Jmol._syncedApplets[i] = applet;
				Jmol._syncedReady[i] = 1;
			} else if (!Jmol._syncedReady[i]) {
				continue;
			}
			n++;
		}
		if (n != Jmol._syncedApplets.length)
			return;
		Jmol._setSyncReady();
	}

	Jmol._setDestroy = function(applet) {
		//MSIE bug responds to any link click even if it is just a JavaScript call

		if (Jmol.featureDetection.allowDestroy)
			Jmol.$windowOn('beforeunload', function () { Jmol._destroy(applet); } );
	}

	Jmol._destroy = function(applet) {
		try {
			if (applet._appletPanel) applet._appletPanel.destroy();
			applet._applet = null;
			Jmol._unsetMouse(applet._canvas)
			applet._canvas = null;
			var n = 0;
			for (var i = 0; i < Jmol._syncedApplets.length; i++) {
				if (Jmol._syncedApplets[i] == applet)
					Jmol._syncedApplets[i] = null;
				if (Jmol._syncedApplets[i])
					n++;
			}
			if (n > 0)
				return;
			Jmol._clearVars();
		} catch(e){}
	}

	////////////// misc core functionality //////////////

	Jmol._setSyncReady = function() {
		Jmol._syncReady = true;
		var s = ""
		for (var i = 0; i < Jmol._syncedApplets.length; i++)
			if (Jmol._syncedCommands[i])
				s += "Jmol.script(Jmol._syncedApplets[" + i + "], Jmol._syncedCommands[" + i + "]);"
		setTimeout(s, 50);  
	}

	Jmol._mySyncCallback = function(appFullName,msg) {
		app = Jmol._applets[appFullName];
		if (app._viewSet) {
			// when can we do this?
//			if (app._viewType == "JSV" && !app._currentView.JMOL)
				Jmol.View.updateFromSync(app, msg);
			return;
		}
		if(!Jmol._syncReady || !Jmol._isJmolJSVSync)
			return 1; // continue processing and ignore me
		for (var i = 0; i < Jmol._syncedApplets.length; i++) {
			if (msg.indexOf(Jmol._syncedApplets[i]._syncKeyword) >= 0)
				Jmol._syncedApplets[i]._syncScript(msg);
		}
		return 0 // prevents further Jmol sync processing 
	}              

	Jmol._getElement = function(applet, what) {
		var d = document.getElementById(applet._id + "_" + what);
		return (d || {});
	} 
	 
	Jmol._evalJSON = function(s,key){
		s = s + "";
		if(!s)
			return [];
		if(s.charAt(0) != "{") {
			if(s.indexOf(" | ") >= 0)
				s = s.replace(/\ \|\ /g, "\n");
			return s;
		}
		var A = (new Function( "return " + s ) )();
		return (!A ? null : key && A[key] != undefined ? A[key] : A);
	}

	Jmol._sortMessages = function(A){
		/*
		 * private function
		 */
		function _sortKey0(a,b){
			return (a[0]<b[0]?1:a[0]>b[0]?-1:0);
		}

		if(!A || typeof (A) != "object")
			return [];
		var B = [];
		for(var i = A.length - 1; i >= 0; i--)
			for(var j = 0, jj= A[i].length; j < jj; j++)
				B[B.length] = A[i][j];
		if(B.length == 0)
			return;
		B = B.sort(_sortKey0);
		return B;
	}

	//////////////////// mouse events //////////////////////

	Jmol._setMouseOwner = function(who, tf) {
		if (who == null || tf) {
			Jmol._mouseOwner = who;
			who && who._applet && (Jmol._lastAppletID = who._applet._id);
		} else if (Jmol._mouseOwner == who) {
			Jmol._mouseOwner = null;
		}
	}

	Jmol._jsGetMouseModifiers = function(ev) {
		var modifiers = 0;
		switch (ev.button) {
		case 0:
			modifiers = 16;//J.api.Event.MOUSE_LEFT;
			break;
		case 1:
			modifiers = 8;//J.api.Event.MOUSE_MIDDLE;
			break;
		case 2:
			modifiers = 4;//J.api.Event.MOUSE_RIGHT;
			break;
		}
		if (ev.shiftKey)
			modifiers += 1;//J.api.Event.SHIFT_MASK;
		if (ev.altKey)
			modifiers += 8;//J.api.Event.ALT_MASK;
		if (ev.ctrlKey)
			modifiers += 2;//J.api.Event.CTRL_MASK;
		return modifiers;
	}

	Jmol._jsGetXY = function(canvas, ev) {
		if (!canvas.applet._ready || Jmol._touching && ev.type.indexOf("touch") < 0)
			return false;
		//ev.preventDefault(); // removed 5/9/2015 -- caused loss of focus on text-box clicking in SwingJS
		var offsets = Jmol.$offset(canvas.id);
		var x, y;
		var oe = ev.originalEvent;
		// drag-drop jQuery event is missing pageX
		ev.pageX || (ev.pageX = oe.pageX);
		ev.pageY || (ev.pageY = oe.pageY);
		Jmol._mousePageX = ev.pageX;
		Jmol._mousePageY = ev.pageY;
		if (oe.targetTouches && oe.targetTouches[0]) {
			x = oe.targetTouches[0].pageX - offsets.left;
			y = oe.targetTouches[0].pageY - offsets.top;
		} else if (oe.changedTouches) {
			x = oe.changedTouches[0].pageX - offsets.left;
			y = oe.changedTouches[0].pageY - offsets.top;
		} else {
			x = ev.pageX - offsets.left;
			y = ev.pageY - offsets.top;
		}
		return (x == undefined ? null : [Math.round(x), Math.round(y), Jmol._jsGetMouseModifiers(ev)]);
	}

  Jmol.setCursor = function(applet, c) {
    if (applet._isJava || applet._readingLocal)
      return;
    var curs;
    switch(c) {
    case 1: 
      curs = "crosshair";
      break;
    case 3: // wait
      curs = "wait";
      Jmol.$setVisible(Jmol.$(applet, "waitimage"), true);
      break;
    case 8: // zoom
      curs = "ns-resize"; 
      break;
    case 12: // hand
      curs = "grab"; 
      break;
    case 13: 
      curs = "move";
      break;
    default:
      Jmol.$setVisible(Jmol.$(applet, "waitimage"), false);
      curs = "default";
      break;
    }
    applet._canvas.style.cursor = curs;
  }
  
	Jmol._gestureUpdate = function(canvas, ev) {
		ev.stopPropagation();
		ev.preventDefault();
		var oe = ev.originalEvent;
		switch (ev.type) {
		case "touchstart":
			Jmol._touching = true;
			break;
		case "touchend":
			Jmol._touching = false;
			break;
		}
		if (!oe.touches || oe.touches.length != 2) return false;
		switch (ev.type) {
		case "touchstart":
			canvas._touches = [[],[]];
			break;
		case "touchmove":
			var offsets = Jmol.$offset(canvas.id);
			var t0 = canvas._touches[0];
			var t1 = canvas._touches[1];
			t0.push([oe.touches[0].pageX - offsets.left, oe.touches[0].pageY - offsets.top]);
			t1.push([oe.touches[1].pageX - offsets.left, oe.touches[1].pageY - offsets.top]);
			var n = t0.length;
			if (n > 3) {
				t0.shift();
				t1.shift();
			}
			if (n >= 2)
				canvas.applet._processGesture(canvas._touches);
			break;
		}
		return true;
	}

	Jmol._jsSetMouse = function(canvas) {

    var doIgnore = function(ev) { return (!ev.target || ("" + ev.target.className).indexOf("swingjs-ui") >= 0) };

		Jmol.$bind(canvas, 'mousedown touchstart', function(ev) {
      if (doIgnore(ev))
        return true;
      if (Jmol._allowKeyboardFocus)	
    	  canvas.focus();
			Jmol._setMouseOwner(canvas, true);
			ev.stopPropagation();
      var ui = ev.target["data-UI"];
      if (!ui || !ui.handleJSEvent(canvas, 501, ev)) 
  			ev.preventDefault();
			canvas.isDragging = true;
			if ((ev.type == "touchstart") && Jmol._gestureUpdate(canvas, ev))
				return !!ui;
			Jmol._setConsoleDiv(canvas.applet._console);
			var xym = Jmol._jsGetXY(canvas, ev);
			if(xym) {
		  	if (ev.button != 2)
          Jmol.Swing.hideMenus(canvas.applet);
        canvas.applet._processEvent(501, xym); //java.awt.Event.MOUSE_DOWN
      }
			return !!ui;
		});
    
		Jmol.$bind(canvas, 'mouseup touchend', function(ev) {
      if (doIgnore(ev))
        return true;
			Jmol._setMouseOwner(null);
			ev.stopPropagation();
      var ui = ev.target["data-UI"];
      if (!ui || !ui.handleJSEvent(canvas, 502, ev))
  			ev.preventDefault();
			canvas.isDragging = false;
			if (ev.type == "touchend" && Jmol._gestureUpdate(canvas, ev))
				return !!ui;
			var xym = Jmol._jsGetXY(canvas, ev);
			if(xym)
  			canvas.applet._processEvent(502, xym);//java.awt.Event.MOUSE_UP
			return !!ui;
		});
    
		Jmol.$bind(canvas, 'mousemove touchmove', function(ev) { // touchmove
      if (doIgnore(ev))
        return true;
      if (Jmol._allowKeyboardFocus)	
    	  canvas.focus(); //in Chrome and Edge, this forces a jump 
      // if the Jmol window is not within viewing bounds
		  // defer to console or menu when dragging within this canvas
			if (Jmol._mouseOwner && Jmol._mouseOwner != canvas && Jmol._mouseOwner.isDragging) {
        if (!Jmol._mouseOwner.mouseMove)
          return true;
	   			Jmol._mouseOwner.mouseMove(ev);
				return false;
			}
			return Jmol._drag(canvas, ev);
		});
		

		Jmol.$bind(canvas, 'keydown keyup', function(ev) {
      if (doIgnore(ev))
        return true;
			ev.stopPropagation();
  			ev.preventDefault();
			var xym = Jmol._jsGetXY(canvas, ev);
			var type = (ev.type == "keydown" ? 401 : 402);
  			canvas.applet._processKeyEvent && canvas.applet._processKeyEvent(type, xym, ev);//java.awt.Event.
			return true;
		});

		Jmol._drag = function(canvas, ev) {
      
			ev.stopPropagation();
			ev.preventDefault();
      
			var isTouch = (ev.type == "touchmove");
			if (isTouch && Jmol._gestureUpdate(canvas, ev))
				return false;
			var xym = Jmol._jsGetXY(canvas, ev);
			if(!xym) return false;
      
			if (!canvas.isDragging)
				xym[2] = 0;

      var ui = ev.target["data-UI"];
      if (canvas.isdragging && (!ui || !ui.handleJSEvent(canvas, 506, ev))) {}
			canvas.applet._processEvent((canvas.isDragging ? 506 : 503), xym); // java.awt.Event.MOUSE_DRAG : java.awt.Event.MOUSE_MOVE
			return !!ui;
}
		
		Jmol.$bind(canvas, 'DOMMouseScroll mousewheel', function(ev) { // Zoom
      if (doIgnore(ev))
        return true;
			ev.stopPropagation();
			ev.preventDefault();
			// Webkit or Firefox
			canvas.isDragging = false;
			var oe = ev.originalEvent;
			var scroll = (oe.detail ? oe.detail : (Jmol.featureDetection.os == "mac" ? 1 : -1) * oe.wheelDelta); // Mac and PC are reverse; but 
			var modifiers = Jmol._jsGetMouseModifiers(ev);
			canvas.applet._processEvent(507,[scroll < 0 ? -1 : 1,0,modifiers]);
			return false;
		});

		// context menu is fired on mouse down, not up, and it's handled already anyway.

		Jmol.$bind(canvas, "contextmenu", function() {return false;});

		Jmol.$bind(canvas, 'mouseout', function(ev) {
      if (doIgnore(ev))
        return true;
      if (Jmol._mouseOwner && !Jmol._mouseOwner.mouseMove) 
        Jmol._setMouseOwner(null);
			if (canvas.applet._appletPanel)
				canvas.applet._appletPanel.startHoverWatcher(false);
			//canvas.isDragging = false;
			var xym = Jmol._jsGetXY(canvas, ev);
			if (!xym)
				return false;
			//canvas.applet._processEvent(502, xym);//J.api.Event.MOUSE_UP
			//canvas.applet._processEvent(505, xym);//J.api.Event.MOUSE_EXITED
			return false;
		});

		Jmol.$bind(canvas, 'mouseenter', function(ev) {
      if (doIgnore(ev))
        return true;
			if (canvas.applet._appletPanel)
				canvas.applet._appletPanel.startHoverWatcher(true);
			if (ev.buttons === 0 || ev.which === 0) {
				canvas.isDragging = false;
				var xym = Jmol._jsGetXY(canvas, ev);
				if (!xym)
					return false;
				canvas.applet._processEvent(504, xym);//J.api.Event.MOUSE_ENTERED	
				canvas.applet._processEvent(502, xym);//J.api.Event.MOUSE_UP
				return false;
			}
		});

	Jmol.$bind(canvas, 'mousemoveoutjsmol', function(evspecial, target, ev) {
      if (doIgnore(ev))
        return true;
		if (canvas == Jmol._mouseOwner && canvas.isDragging) {
			return Jmol._drag(canvas, ev);
		}
	});

		if (canvas.applet._is2D)
			Jmol.$resize(function() {
				if (!canvas.applet)
					return;
				canvas.applet._resize();
			});
 
		Jmol.$bind('body', 'mouseup touchend', function(ev) {
      if (doIgnore(ev))
        return true;
			if (canvas.applet)
				canvas.isDragging = false;
			Jmol._setMouseOwner(null);
		});

	}

	Jmol._jsUnsetMouse = function(canvas) {
		canvas.applet = null;
		Jmol.$bind(canvas, 'mousedown touchstart mousemove touchmove mouseup touchend DOMMouseScroll mousewheel contextmenu mouseout mouseenter', null);
		Jmol._setMouseOwner(null);
	}


////// Jmol.Swing interface  for Javascript implementation of Swing dialogs and menus

Jmol.Swing = {
	// a static class
	count:0,
	menuInitialized:0,
	menuCounter:0,
	htDialogs:{}
};

(function(Swing) {

SwingController = Swing;

Swing.setDraggable = function(Obj) {
	
	var proto = Obj.prototype;
	if (proto.setContainer)
		return;
	
	// for menus, console, and 
	proto.setContainer = function(container) {
		this.container = container;
		container.obj = this;
		this.isDragging = false;
		this.ignoreMouse = false;
		var me = this;
		container.bind('mousedown touchstart', function(ev) {
			if (me.ignoreMouse) {
				me.ignoreMouse = false;
				return true;
			}
			Jmol._setMouseOwner(me, true);
			me.isDragging = true;
			me.pageX = ev.pageX;
			me.pageY = ev.pageY;
			return false;
		});
		container.bind('mousemove touchmove', function(ev) {
			if (me.isDragging && Jmol._mouseOwner == me) {
				me.mouseMove(ev);
				return false;
			}
		});
		container.bind('mouseup touchend', function(ev) {
			me.mouseUp(ev);
			Jmol._setMouseOwner(null);
		});
	};

	proto.mouseUp = function(ev) {
		if (this.isDragging && Jmol._mouseOwner == this) {
			this.pageX0 += (ev.pageX - this.pageX);
			this.pageY0 += (ev.pageY - this.pageY);
			this.isDragging = false;
			return false;
		}
		Jmol._setMouseOwner(null);
	}

	proto.setPosition = function() {
		if (Jmol._mousePageX === null) {
			var id = this.applet._id + "_" + (this.applet._is2D ? "canvas2d" : "canvas");
			var offsets = Jmol.$offset(id);
			Jmol._mousePageX = offsets.left;
			Jmol._mousePageY = offsets.top;
		}
		this.pageX0 = Jmol._mousePageX;
		this.pageY0 = Jmol._mousePageY;
		var pos = { top: Jmol._mousePageY + 'px', left: Jmol._mousePageX + 'px' };
		this.container.css(pos);
	};

	proto.mouseMove = function(ev) {
		if (this.isDragging && Jmol._mouseOwner == this) {
			this.timestamp = System.currentTimeMillis(); // used for menu closure
			var x = this.pageX0 + (ev.pageX - this.pageX);
			var y = this.pageY0 + (ev.pageY - this.pageY);
      Jmol._mousePageX = x;
      Jmol._mousePageY = y;
			this.container.css({ top: y + 'px', left: x + 'px' })
		}
	};

	proto.dragBind = function(isBind) {
		this.applet._ignoreMouse = !isBind;
		this.container.unbind('mousemoveoutjsmol');
		this.container.unbind('touchmoveoutjsmol');
		this.container.unbind('mouseupoutjsmol');
		this.container.unbind('touchendoutjsmol');
		Jmol._setMouseOwner(null);
		if (isBind) {
			var me = this;
			this.container.bind('mousemoveoutjsmol touchmoveoutjsmol', function(evspecial, target, ev) {
				me.mouseMove(ev);
			});
			this.container.bind('mouseupoutjsmol touchendoutjsmol', function(evspecial, target, ev) {
				me.mouseUp(ev);
			});
		}
	};
}

// Dialog //

Swing.JSDialog = function () {
}

Swing.setDraggable(Swing.JSDialog);

///// calls from javajs and other Java-derived packages /////

Swing.getScreenDimensions = function(d) {
	d.width = $(window).width();
	d.height = $(window).height();
}

Swing.dispose = function(dialog) {
	Jmol.$remove(dialog.id + "_mover");
	delete Swing.htDialogs[dialog.id]
	dialog.container.obj.dragBind(false);
//  var btns = $("#" + dialog.id + " *[id^='J']"); // add descendents with id starting with "J"
//  for (var i = btns.length; --i >= 0;)
//    delete Dialog.htDialogs[btns[i].id]
	//System.out.println("JSmolCore.js: dispose " + dialog.id)
}
 
Swing.register = function(dialog, type) {
	dialog.id = type + (++Swing.count);
	Swing.htDialogs[dialog.id] = dialog;
	//System.out.println("JSmolCore.js: register " + dialog.id)

}

Swing.setDialog = function(dialog) {
	Jmol._setMouseOwner(null);
	Jmol.$remove(dialog.id);
	//System.out.println("removed " + dialog.id)
	var id = dialog.id + "_mover";
	var container = Jmol._$(id);
	var jd;
	//System.out.println("JSmolCore.js: setDialog " + dialog.id);
	if (container[0]) {
		container.html(dialog.html);
		jd = container[0].jd;
	} else {
		Jmol.$after("body","<div id='" + id + "' style='position:absolute;left:0px;top:0px;'>" + dialog.html + "</div>");
		var jd = new Swing.JSDialog();
		container = Jmol._$(id);
		dialog.container = container;
		jd.applet = dialog.manager.vwr.html5Applet;
		jd.setContainer(container);
		jd.dialog = dialog;
		jd.setPosition();  
		jd.dragBind(true);
		container[0].jd = jd; 
	}
	Jmol.$bind("#" + dialog.id + " .JButton", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JComboBox", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JCheckBox", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JTextField", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JTable", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JScrollPane", "mousedown touchstart", function(event) { jd.ignoreMouse=true });
	Jmol.$bind("#" + dialog.id + " .JEditorPane", "mousedown touchstart", function(event) { jd.ignoreMouse=true });

}
 
Swing.setSelected = function(chk) {
 Jmol.$prop(chk.id, 'checked', !!chk.selected);
}

Swing.setSelectedIndex = function(cmb) {
 Jmol.$prop(cmb.id, 'selectedIndex', cmb.selectedIndex);
}

Swing.setText = function(btn) {
 Jmol.$prop(btn.id, 'value', btn.text);
}

Swing.setVisible = function(c) {
	Jmol.$setVisible(c.id, c._visible);
}

Swing.setEnabled = function(c) {
	Jmol.$setEnabled(c.id, c.enabled);
}

/// callbacks from the HTML elements ////
 
Swing.click = function(element, keyEvent) {
	var component = Swing.htDialogs[element.id];
	if (component) {
		//System.out.println("click " + element + " " + component)
		var info = component.toString();
		// table cells will have an id but are not registered
		if (info.indexOf("JCheck") >= 0) {
				component.selected = element.checked;
		} else if (info.indexOf("JCombo") >= 0) {
			component.selectedIndex = element.selectedIndex;
		} else if (component.text != null) {  // JButton, JTextField
			component.text = element.value;
			if (keyEvent && ((keyEvent.charCode || keyEvent.keyCode) != 13))
				return;
		}    
	}
	var dialog = Swing.htDialogs[Jmol.$getAncestorDiv(element.id, "JDialog").id];
	var key = (component ? component.name :  dialog.registryKey + "/" + element.id);
	//System.out.println("JSmolCore.js: click " + key); 
	dialog.manager.actionPerformed(key);
}

Swing.setFront = function(dialog) {
  var applet = dialog.manager.vwr.html5Applet;
	if (dialog.zIndex != Jmol._getZ(applet, "dialog"))
	 dialog.zIndex = Jmol._incrZ(applet, "dialog");
	dialog.container && ((dialog.container[0] || dialog.container).style.zIndex = dialog.zIndex);
}

Swing.hideMenus = function(applet) {
	// called from LEFT-DOWN mouse event
	var menus = applet._menus;
	if (menus)
		for (var i in menus)
			if (menus[i]._visible)
				Swing.hideMenu(menus[i]);
}

Swing.windowClosing = function(element) {
	var dialog = Swing.htDialogs[Jmol.$getAncestorDiv(element.id, "JDialog").id];
	if (dialog.registryKey) {
		//System.out.println("JSmolCore.js: windowClosing " + dialog.registryKey); 
		dialog.manager.processWindowClosing(dialog.registryKey);
	} else {
		//System.out.println("JSmolCore.js: windowClosing " + dialog.title); 
		dialog.dispose();
	}
}

})(Jmol.Swing);

Jmol._track = function(applet) {
	// this function inserts an iFrame that can be used to track your page's applet use. 
	// By default it tracks to a page at St. Olaf College, but you can change that. 
	// and you can use
	//
	// delete Jmol._tracker
	//
	// yourself to not have you page execute this 
	//
	if (Jmol._tracker){
		try {  
			var url = Jmol._tracker + "&applet=" + applet._jmolType + "&version=" + Jmol._version 
				+ "&appver=" + Jmol.___JmolVersion + "&url=" + encodeURIComponent(document.location.href);
			var s = '<iframe style="display:none" width="0" height="0" frameborder="0" tabindex="-1" src="' + url + '"></iframe>'
			Jmol.$after("body", s);
		} catch (e) {
			// ignore
		}
		delete Jmol._tracker;
	}
	return applet;
}

var __profiling;

Jmol.getProfile = function(doProfile) {
  if (!self.Clazz || !self.JSON)
    return;
  if (!__profiling)
    Clazz._startProfiling(__profiling = (arguments.length == 0 || doProfile));
	return Clazz.getProfile();
}

Jmol._getInChIKey = function(applet, data) {
	if (data.indexOf("MOL=") >= 0)
		data = data.split("MOL=")[1].split("\"")[0];

}

Jmol._getAttr = function(s, a) {
	var pt = s.indexOf(a + "=");
	return (pt >= 0 && (pt = s.indexOf('"', pt)) >= 0 
		? s.substring(pt+1, s.indexOf('"', pt+1)) : null);
}

Jmol.User = {
	viewUpdatedCallback: null
}

Jmol.View = {

// The objective of Jmol.View is to coordinate
// asynchronous applet loading and atom/peak picking
// among any combination of Jmol, JME, and JSV.
// 
// basic element is a view object:
//   view = {
//     viewType1: viewRecord1,
//     viewType2: viewRecord2,
//     viewType3: viewRecord3
//   }
// where viewType is one of (Jmol, JME, JSV)
// and a viewRecord is an object
// with elements .chemID, .applet, .data, .script
//
// Jmol.View.sets is a list of cached views[0..n]
// for a given group of applet objects with common Info.viewSet
//
// Bob Hanson 1/22/2014 7:05:38 AM

	count: 0,
	applets: {},
	sets: {}
};

(function(View) {

// methods called from other modules have no "_" in their name

View.resetView = function(applet, appletNot) {
  debugger;
  if (appletNot) {
  	if (!appletNot._viewSet)
  		return;
    var set = Jmol.View.applets[appletNot._viewSet]
    for (var applet in set) {
      if (applet == appletNot)
        continue;
      Jmol.View.resetView(applet);
    }
    return;
  }
  if (applet) {
  	applet._reset();
    Jmol.View.updateView(applet);
  }
}

View.updateView = function(applet, Info, _View_updateView) {
	// Info.chemID, Info.data, possibly Info.viewID if no chemID
	// return from applet after asynchronous download of new data
	if (applet._viewSet == null)
		return;
  Info || (Info = {});
	Info.chemID || (applet._searchQuery = null);
	Info.data || (Info.data = "N/A");
	Info.type = applet._viewType;
	if((applet._currentView = View.__findView(applet._viewSet, Info)) == null)
		applet._currentView = View.__createViewSet(applet._viewSet, Info.chemID, Info.viewID || Info.chemID);
	applet._currentView[Info.type].data = Info.data;
	applet._currentView[Info.type].smiles = applet._getSmiles();
	if (Jmol.User.viewUpdatedCallback)
		Jmol.User.viewUpdatedCallback(applet, "updateView");
	View.__setView(applet._currentView, applet, false);
}

View.updateFromSync = function(applet, msg) {
	applet._updateMsg = msg;
	var id = Jmol._getAttr(msg, "sourceID") || Jmol._getAttr(msg, "file");
	if (!id)
		return;
	var view = View.__findView(applet._viewSet, {viewID:id});
	if (view == null)
		return Jmol.updateView(applet, msg); // JSV has been updated internally
	if (view != applet._currentView)
		View.__setView(view, applet, true);
	var A = ((id = Jmol._getAttr(msg, "atoms")) && msg.indexOf("selectionhalos ON") >= 0  
		? eval("[" + id + "]") : []);
	setTimeout(function(){if (applet._currentView == view)View.updateAtomPick(applet, A)}, 10);
	if (Jmol.User.viewUpdatedCallback)
		Jmol.User.viewUpdatedCallback(applet, "updateFromSync");
}

View.updateAtomPick = function(applet, A, _View_updateAtomPick) {
	var view = applet._currentView;
	if (view == null)
		return;
	for (var viewType in view) {
		if (viewType != "info" && view[viewType].applet != applet)
			view[viewType].applet._updateAtomPick(A);
	}
	if (Jmol.User.viewUpdatedCallback)
		Jmol.User.viewUpdatedCallback(applet, "updateAtomPick");
}

View.dumpViews = function(setID) {
	var views = View.sets[setID];
	if (!views)
	  return;
	var s = "View set " + setID + ":\n";
	var applets = View.applets[setID];
	for (var i in applets)
		s += "\napplet " + applets[i]._id
			+ " currentView=" + (applets[i]._currentView ? applets[i]._currentView.info.viewID : null);
	for (var i = views.length; --i >= 0;) {
		var view = views[i];
		s += "\n\n<b>view=" + i 
			+ " viewID=" + view.info.viewID 
			+ " chemID=" + view.info.chemID + "</b>\n"
		var v;
		for (var viewType in view) 
			if (viewType != "info")
				s += "\nview=" + i + " type=" + viewType + " applet=" 
					+ ((v = view[viewType]).applet ? v.applet._id : null) 
					+ " SMILES=" + v.smiles + "\n"
					+ " atomMap=" + JSON.stringify(v.atomMap) + "\n"
					+ " data=\n" + v.data + "\n"
	}
	return s
}


// methods starting with "__" are "private" to JSmolCore.js

View.__init = function(applet) {
  var set = applet._viewSet;
	var a = View.applets;
	a[set] || (a[set] = {});
	a[set][applet._viewType] = applet;
}

View.__findView = function(set, Info) {
	var views = View.sets[set];
	if (views == null)
		views = View.sets[set] = [];
	for (var i = views.length; --i >= 0;) {
		var view = views[i];
		if (Info.viewID) {
			if (view.info.viewID == Info.viewID)
				return view;
		} else if (Info.chemID != null && Info.chemID == view.info.chemID) {
				return view;
		} else {
			for (var viewType in view) { 
				if (viewType != "info") {
					if (Info.data != null && view[viewType].data != null ? Info.data == view[viewType].data 
						: Info.type == viewType)
							return view;
				}
			}
		}
	}
	return null;  
}

View.__createViewSet = function(set, chemID, viewID, _createViewSet) {
	View.count++;
	var view = {info:{chemID: chemID, viewID: viewID || "model_" + View.count}};

	for (var id in Jmol._applets) {
		var a = Jmol._applets[id];
		if (a._viewSet == set)
			view[a._viewType] = {applet:a, data: null};
	}
	View.sets[set].push(view);
	return view;
}

View.__setView = function(view, applet, isSwitch, _setView) {
	// called from Jmol._searchMol and Jmol.View.setCurrentView 
	// to notify the applets in the set that there may be new data for them
	// skip the originating applet itself and cases where the data has not changed.
	// stop at first null data, because that will initiate some sort of asynchronous
	// call that will be back here afterward.

	for (var viewType in view) {
			if (viewType == "info") 
				continue;
		var rec = view[viewType];
		var a = rec.applet;
		var modified = (isSwitch || a != null && a._molData == "<modified>");

		if (a == null || a == applet && !modified)
			continue; // may be a mol3d required by JSV but not having a corresponding applet
		var wasNull = (rec.data == null);
		var haveView = (a._currentView != null);
		a._currentView = view; 
		if (haveView && view[viewType].data == rec.data && !wasNull & !modified)
			continue;
		a._loadModelFromView(view);
		if (wasNull)
			break;
	}

	// Either all are taken care of or one was null,
	// in which case we have started an asynchronous
	// process to get the data, and we can quit here.
	// In either case, we are done.
}

}) (Jmol.View);

Jmol.Cache = {fileCache: {}};

Jmol.Cache.get = function(filename) {
	return Jmol.Cache.fileCache[filename];
}

Jmol.Cache.put = function(filename, data) {
  Jmol.Cache.fileCache[filename] = data;
}

	Jmol.Cache.setDragDrop = function(me, divname) {
		Jmol.$appEvent(me, divname, "dragover", function(e) {
			e = e.originalEvent;
			e.stopPropagation();
			e.preventDefault();
			e.dataTransfer.dropEffect = 'copy';
		});
		Jmol.$appEvent(me, divname, "drop", function(e) {
			var oe = e.originalEvent;
			oe.stopPropagation();
			oe.preventDefault();
			var file = oe.dataTransfer.files[0];
			if (file == null) {
				// FF and Chrome will drop an image here
				// but it will be only a URL, not an actual file. 
				try {
				  file = "" + oe.dataTransfer.getData("text");
				  if (file.indexOf("file:/") == 0 || file.indexOf("http:/") == 0 || file.indexOf("https:/") == 0) {
				  	me._scriptLoad(file);
				  	return;
			  	}
				} catch(e) {
				  return;
				}
			  // some other format
			  return;
			}
			// MSIE will drop an image this way, though, and load it!
			var reader = new FileReader();
			reader.onloadend = function(evt) {
				if (evt.target.readyState == FileReader.DONE) {
					var cacheName = "cache://DROP_" + file.name;
					var bytes = Clazz.newByteArray(-1, evt.target.result);
					if (divname == "console_input") {
						var s = String.instantialize(bytes);
						Jmol.$(me,divname).val(s.indexOf('\0') < 0 ? s : bytes.length + " bytes");
						return;
					}
					if (!cacheName.endsWith(".spt"))
						me._appletPanel.cacheFileByName("cache://DROP_*",false);
					if (me._viewType == "JSV" || cacheName.endsWith(".jdx")) // shared by Jmol and JSV
						Jmol.Cache.put(cacheName, bytes);
					else
						me._appletPanel.cachePut(cacheName, bytes);
					var xym = Jmol._jsGetXY(me._canvas, e);
					var retType = [null];
					if(xym && (!me._appletPanel.setStatusDragDropped || me._appletPanel.setStatusDragDropped(0, xym[0], xym[1], cacheName, retType))) {
						me._appletPanel.openFileAsyncSpecialType(cacheName, 1, retType[0]);
					}
				}
			};
			reader.readAsArrayBuffer(file);
		});
	}
  
})(Jmol, jQuery);
