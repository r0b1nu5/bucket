window.html10n=function(window,document,undefined){function Loader(n){this.resources=n,this.cache={},this.langs={}}function getPluralRules(n){function t(n,t){return-1!==t.indexOf(n)}function e(n,t,e){return n>=t&&e>=n}var r={af:3,ak:4,am:4,ar:1,asa:3,az:0,be:11,bem:3,bez:3,bg:3,bh:4,bm:0,bn:3,bo:0,br:20,brx:3,bs:11,ca:3,cgg:3,chr:3,cs:12,cy:17,da:3,de:3,dv:3,dz:0,ee:3,el:3,en:3,eo:3,es:3,et:3,eu:3,fa:0,ff:5,fi:3,fil:4,fo:3,fr:5,fur:3,fy:3,ga:8,gd:24,gl:3,gsw:3,gu:3,guw:4,gv:23,ha:3,haw:3,he:2,hi:4,hr:11,hu:0,id:0,ig:0,ii:0,is:3,it:3,iu:7,ja:0,jmc:3,jv:0,ka:0,kab:5,kaj:3,kcg:3,kde:0,kea:0,kk:3,kl:3,km:0,kn:0,ko:0,ksb:3,ksh:21,ku:3,kw:7,lag:18,lb:3,lg:3,ln:4,lo:0,lt:10,lv:6,mas:3,mg:4,mk:16,ml:3,mn:3,mo:9,mr:3,ms:0,mt:15,my:0,nah:3,naq:7,nb:3,nd:3,ne:3,nl:3,nn:3,no:3,nr:3,nso:4,ny:3,nyn:3,om:3,or:3,pa:3,pap:3,pl:13,ps:3,pt:3,rm:3,ro:9,rof:3,ru:11,rwk:3,sah:0,saq:3,se:7,seh:3,ses:0,sg:0,sh:11,shi:19,sk:12,sl:14,sma:7,smi:7,smj:7,smn:7,sms:7,sn:3,so:3,sq:3,sr:11,ss:3,ssy:3,st:3,sv:3,sw:3,syr:3,ta:3,te:3,teo:3,th:0,ti:4,tig:3,tk:3,tl:4,tn:3,to:0,tr:0,ts:3,tzm:22,uk:11,ur:3,ve:3,vi:0,vun:3,wa:4,wae:3,wo:0,xh:3,xog:3,yo:0,zh:0,zu:3},o={0:function(n){return"other"},1:function(n){return e(n%100,3,10)?"few":0===n?"zero":e(n%100,11,99)?"many":2==n?"two":1==n?"one":"other"},2:function(n){return 0!==n&&n%10===0?"many":2==n?"two":1==n?"one":"other"},3:function(n){return 1==n?"one":"other"},4:function(n){return e(n,0,1)?"one":"other"},5:function(n){return e(n,0,2)&&2!=n?"one":"other"},6:function(n){return 0===n?"zero":n%10==1&&n%100!=11?"one":"other"},7:function(n){return 2==n?"two":1==n?"one":"other"},8:function(n){return e(n,3,6)?"few":e(n,7,10)?"many":2==n?"two":1==n?"one":"other"},9:function(n){return 0===n||1!=n&&e(n%100,1,19)?"few":1==n?"one":"other"},10:function(n){return e(n%10,2,9)&&!e(n%100,11,19)?"few":n%10!=1||e(n%100,11,19)?"other":"one"},11:function(n){return e(n%10,2,4)&&!e(n%100,12,14)?"few":n%10===0||e(n%10,5,9)||e(n%100,11,14)?"many":n%10==1&&n%100!=11?"one":"other"},12:function(n){return e(n,2,4)?"few":1==n?"one":"other"},13:function(n){return e(n%10,2,4)&&!e(n%100,12,14)?"few":1!=n&&e(n%10,0,1)||e(n%10,5,9)||e(n%100,12,14)?"many":1==n?"one":"other"},14:function(n){return e(n%100,3,4)?"few":n%100==2?"two":n%100==1?"one":"other"},15:function(n){return 0===n||e(n%100,2,10)?"few":e(n%100,11,19)?"many":1==n?"one":"other"},16:function(n){return n%10==1&&11!=n?"one":"other"},17:function(n){return 3==n?"few":0===n?"zero":6==n?"many":2==n?"two":1==n?"one":"other"},18:function(n){return 0===n?"zero":e(n,0,2)&&0!==n&&2!=n?"one":"other"},19:function(n){return e(n,2,10)?"few":e(n,0,1)?"one":"other"},20:function(n){return!e(n%10,3,4)&&n%10!=9||e(n%100,10,19)||e(n%100,70,79)||e(n%100,90,99)?n%1e6===0&&0!==n?"many":n%10!=2||t(n%100,[12,72,92])?n%10!=1||t(n%100,[11,71,91])?"other":"one":"two":"few"},21:function(n){return 0===n?"zero":1==n?"one":"other"},22:function(n){return e(n,0,1)||e(n,11,99)?"one":"other"},23:function(n){return e(n%10,1,2)||n%20===0?"one":"other"},24:function(n){return e(n,3,10)||e(n,13,19)?"few":t(n,[2,12])?"two":t(n,[1,11])?"one":"other"}},i=r[n.replace(/-.*$/,"")];return i in o?o[i]:(console.warn("plural form unknown for ["+n+"]"),function(){return"other"})}function asyncForEach(n,t,e){var r=0,o=n.length;t(n[r],r,function i(a){return a&&console.error(a),r++,o>r?t(n[r],r,i):void e()})}function getTranslatableChildren(n){if(!document.querySelectorAll){if(!n)return[];for(var t=n.getElementsByTagName("*"),e=[],r=0,o=t.length;o>r;r++)t[r].getAttribute("data-l10n-id")&&e.push(t[r]);return e}return n.querySelectorAll("*[data-l10n-id]")}function substArguments(n,t){for(var e,r=/\{\{\s*([a-zA-Z\.]+)\s*\}\}/;e=r.exec(n);){if(!e||e.length<2)return n;var o=e[1],i="";if(o in t)i=t[o];else{if(!(o in translations))return console.warn("Could not find argument {{"+o+"}}"),n;i=translations[o]}n=n.substring(0,e.index)+i+n.substr(e.index+e[0].length)}return n}function substMacros(n,t,e){for(var r,o=/\{\[\s*([a-zA-Z]+)\(([a-zA-Z]+)\)((\s*([a-zA-Z]+)\: ?([ a-zA-Z{}]+),?)+)*\s*\]\}/;r=o.exec(t);){var i=r[1],a=r[2],s=r[3],u={};if(i in html10n.macros){s&&s.match(/(?=\s*)([a-zA-Z]+)\: ?([ a-zA-Z{}]+)(?=,?)/g).forEach(function(n){var t=n.split(":"),e=t[0],r=t[1].trim();u[e]=r});var l;e&&a in e?l=e[a]:a in html10n.translations&&(l=translations[a]);var f=html10n.macros[i];t=t.substr(0,r.index)+f(n,l,u)+t.substr(r.index+r[0].length)}}return t}!function(){for(var n=function(){},t=["log","debug","info","warn","error","assert","dir","dirxml","group","groupEnd","time","timeEnd","count","trace","profile","profileEnd"],e=window.console=window.console||{},r=0;r<t.length;++r)e[t[r]]||(e[t[r]]=n)}(),Array.prototype.forEach||(Array.prototype.forEach=function(n,t){for(var e=0,r=this.length;r>e;++e)e in this&&n.call(t,this[e],e,this)}),Array.prototype.indexOf||(Array.prototype.indexOf=function(n){"use strict";if(null==this)throw new TypeError;var t=Object(this),e=t.length>>>0;if(0===e)return-1;var r=0;if(arguments.length>1&&(r=Number(arguments[1]),r!=r?r=0:0!=r&&r!=1/0&&r!=-(1/0)&&(r=(r>0||-1)*Math.floor(Math.abs(r)))),r>=e)return-1;for(var o=r>=0?r:Math.max(e-Math.abs(r),0);e>o;o++)if(o in t&&t[o]===n)return o;return-1});var MicroEvent=function(){};MicroEvent.prototype={bind:function(n,t){this._events=this._events||{},this._events[n]=this._events[n]||[],this._events[n].push(t)},unbind:function(n,t){this._events=this._events||{},n in this._events!=!1&&this._events[n].splice(this._events[n].indexOf(t),1)},trigger:function(n){if(this._events=this._events||{},n in this._events!=!1)for(var t=0;t<this._events[n].length;t++)this._events[n][t].apply(this,Array.prototype.slice.call(arguments,1))}},MicroEvent.mixin=function(n){var t=["bind","unbind","trigger"];if(n)for(var e=0;e<t.length;e++)n[t[e]]=MicroEvent.prototype[t[e]]},Loader.prototype.load=function(n,t){if(this.langs[n])return t();if(this.resources.length>0)for(var e=0,r=0,o=this.resources.length;o>r;r++)this.fetch(this.resources[r],n,function(n){e++,n&&console.warn(n),o>e||t&&t()})},Loader.prototype.fetch=function(n,t,e){var r=this;if(this.cache[n])return void this.parse(t,n,this.cache[n],e);var o=new XMLHttpRequest;o.open("GET",n,!0),o.overrideMimeType&&o.overrideMimeType("application/json; charset=utf-8"),o.onreadystatechange=function(){if(4==o.readyState)if(200==o.status||0===o.status){var i=JSON.parse(o.responseText);r.cache[n]=i,r.parse(t,n,i,e)}else e(new Error("Failed to load "+n))},o.send(null)},Loader.prototype.parse=function(n,t,e,r){if("object"!=typeof e)return void r(new Error("A file couldn't be parsed as json."));if(!e[n]){var o,i="Couldn't find translations for "+n;~n.indexOf("-")&&(n=n.split("-")[0]);for(o in e)if(n!=o&&0===o.indexOf(n)&&e[o]){n=o;break}if(n!=o)return r(new Error(i))}if("string"==typeof e[n]){var a=e[n];return 0!=e[n].indexOf("http")&&0!=e[n].indexOf("/")&&(a=t+"/../"+e[n]),void this.fetch(a,n,r)}return"object"!=typeof e[n]?void r(new Error("Translations should be specified as JSON objects!")):(this.langs[n]=e[n],void r())};var html10n={language:null};return MicroEvent.mixin(html10n),html10n.macros={},html10n.rtl=["ar","dv","fa","ha","he","ks","ku","ps","ur","yi"],html10n.macros.plural=function(n,t,e){var r,o=parseFloat(t);if(!isNaN(o)){this._pluralRules||(this._pluralRules=getPluralRules(html10n.language));var i=this._pluralRules(o);return 0===o&&"zero"in e?r=e.zero:1==o&&"one"in e?r=e.one:2==o&&"two"in e?r=e.two:i in e&&(r=e[i]),r}},html10n.localize=function(n){var t=this;"string"==typeof n&&(n=[n]);var e=0;n.forEach(function(t){t&&(n[e++]=t,~t.indexOf("-")&&(n[e++]=t.substr(0,t.indexOf("-"))))}),this.build(n,function(n,e){html10n.translations=e,html10n.translateElement(e),t.trigger("localized")})},html10n.translateElement=function(n,t){t=t||document.documentElement;for(var e=t?getTranslatableChildren(t):document.childNodes,r=0,o=e.length;o>r;r++)this.translateNode(n,e[r]);this.translateNode(n,t)},html10n.get=function(n,t){var e=html10n.translations;if(!e)return console.warn("No translations available (yet)");if(!e[n])return console.warn("Could not find string "+n);var r=e[n];return r=substMacros(n,r,t),r=substArguments(r,t)},html10n.translateNode=function(translations,node){var str={};if(str.id=node.getAttribute("data-l10n-id"),str.id){if(!translations[str.id])return console.warn("Couldn't find translation key "+str.id);if(window.JSON)str.args=JSON.parse(node.getAttribute("data-l10n-args"));else try{str.args=eval(node.getAttribute("data-l10n-args"))}catch(e){console.warn("Couldn't parse args for "+str.id)}str.str=html10n.get(str.id,str.args);var prop,index=str.id.lastIndexOf("."),attrList={title:1,innerHTML:1,alt:1,textContent:1,value:1,placeholder:1};if(prop=index>0&&str.id.substr(index+1)in attrList?str.id.substr(index+1):document.body.textContent?"textContent":"innerText",0===node.children.length||"textContent"!=prop)node[prop]=str.str,node.setAttribute("aria-label",str.str);else{for(var children=node.childNodes,found=!1,i=0,n=children.length;n>i;i++)3===children[i].nodeType&&/\S/.test(children[i].textContent)&&(found?children[i].nodeValue="":(children[i].nodeValue=str.str,found=!0));found||console.warn("Unexpected error: could not translate element content for key "+str.id,node)}}},html10n.build=function(n,t){var e=this,r={};asyncForEach(n,function(n,t,r){return n?void e.loader.load(n,r):r()},function(){var o;n.reverse();for(var i=0,a=n.length;a>i;i++)if(o=n[i]){if(!(o in e.loader.langs)){~o.indexOf("-")&&(o=o.split("-")[0]);for(var s in e.loader.langs)if(o!=s&&0===s.indexOf(o)){o=s;break}if(o!=s)continue}for(var u in e.loader.langs[o])r[u]=e.loader.langs[o][u];e.language=o}t(null,r)})},html10n.getLanguage=function(){return this.language},html10n.getDirection=function(){if(this.language){var n=-1==this.language.indexOf("-")?this.language:this.language.substr(0,this.language.indexOf("-"));return-1==html10n.rtl.indexOf(n)?"ltr":"rtl"}},html10n.index=function(){for(var n=document.getElementsByTagName("link"),t=[],e=0,r=n.length;r>e;e++)"application/l10n+json"==n[e].type&&t.push(n[e].href);this.loader=new Loader(t),this.trigger("indexed")},document.addEventListener?document.addEventListener("DOMContentLoaded",function(){html10n.index()},!1):window.attachEvent&&window.attachEvent("onload",function(){html10n.index()},!1),window._===undefined&&(window._=html10n.get),html10n}(window,document);