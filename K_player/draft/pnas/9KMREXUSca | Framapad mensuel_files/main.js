$(document).ready(function(){var e=new RegExp(/.*\/p\/[^\/]+/).exec(document.location.pathname)||clientVars.padId;document.location.href.replace(document.location.pathname,e);$("#exportmarkdowna").attr("href",e+"/export/markdown")});