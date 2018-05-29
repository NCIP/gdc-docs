/******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__webpack_require__.d = function(exports, name, getter) {
/******/ 		if(!__webpack_require__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, {
/******/ 				configurable: false,
/******/ 				enumerable: true,
/******/ 				get: getter
/******/ 			});
/******/ 		}
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__webpack_require__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__webpack_require__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(__webpack_require__.s = 3);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
var baseUrl = 'https://gdc-mvs.nci.nih.gov/gdc/search';

var api = {
  suggest: function suggest(value, callback) {
    $.getJSON({
      url: baseUrl + "/suggest?keyword=" + value,
      success: function success(data) {
        callback(data);
      }
    });
  },
  searchAll: function searchAll(keyword, option, callback) {
    $.getJSON(baseUrl + '/all/p', { keyword: keyword, option: JSON.stringify(option) }, function (result) {
      var items = result;
      callback(keyword, option, items);
    });
  },
  getGDCDataById: function getGDCDataById(id, callback) {
    $.getJSON(baseUrl + '/p/local/vs', { id: id }, function (result) {
      callback(id, result);
    });
  },
  getCDEDataById: function getCDEDataById(id, callback) {
    $.getJSON(baseUrl + '/p/cde/vs', { id: id }, function (result) {
      callback(id, result);
    });
  },
  getGDCandCDEDataById: function getGDCandCDEDataById(ids, callback) {
    $.getJSON(baseUrl + '/p/both/vs', { local: ids.local, cde: ids.cde }, function (result) {
      callback(ids, result);
    });
  },
  evsRestApi: function evsRestApi(id, callback) {
    $.getJSON(baseUrl + '/ncit/detail?code=' + id, function (result) {
      callback(id, result);
    });
  }
};

/* harmony default export */ __webpack_exports__["a"] = (api);

/***/ }),
/* 1 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

var _prevScrollOffset = 0;
var heightSlider = $('.navbar .container').height();
var windowEl = $(window);
var _headerOffset = heightSlider;

function _onScroll() {
  var currentScrollOffset = windowEl.scrollTop();
  var delta = currentScrollOffset - _prevScrollOffset;

  if (delta > 0) {
    _headerOffset = heightSlider - 64;
    _prevScrollOffset = currentScrollOffset;
  } else {
    _headerOffset = heightSlider;
    _prevScrollOffset = currentScrollOffset;
  }
}

function _onResize() {
  heightSlider = $('.navbar .container').height();
  _headerOffset = heightSlider;
}

windowEl.scroll(_onScroll);

windowEl.resize(_onResize);

var func = {
  headerOffset: function headerOffset() {
    return _headerOffset;
  }
};

/* harmony default export */ __webpack_exports__["a"] = (func);

/***/ }),
/* 2 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony export (immutable) */ __webpack_exports__["a"] = render;
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__result_table___ = __webpack_require__(6);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__props_table___ = __webpack_require__(8);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__values_table___ = __webpack_require__(10);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_3__tabs___ = __webpack_require__(12);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_4__shared__ = __webpack_require__(1);






function render(keyword, option, items) {
  var html = "";
  if (items.length !== 0) {
    var trsHtml = __WEBPACK_IMPORTED_MODULE_0__result_table___["a" /* default */].render(items);
    trsHtml.active = false;
    var psHtml = __WEBPACK_IMPORTED_MODULE_1__props_table___["a" /* default */].render(items);
    psHtml.active = false;
    var vsHtml = __WEBPACK_IMPORTED_MODULE_2__values_table___["a" /* default */].render(items);
    vsHtml.active = false;
    if (option.activeTab == 0) {
      vsHtml.active = true;
    } else if (option.activeTab == 1) {
      psHtml.active = true;
    } else {
      trsHtml.active = true;
    }
    html = Object(__WEBPACK_IMPORTED_MODULE_3__tabs___["a" /* default */])(trsHtml, psHtml, vsHtml, keyword);
  } else if (option.error == true) {
    html = '<div class="indicator indicator--has-error">Please, enter a valid keyboard!</div>';
  } else {
    html = '<div class="indicator">Sorry, no results found for kerword: <span class="indicator__term">' + keyword + '</span></div>';
  }

  $("#root").html(html);

  if ($("#tree_table").length) {
    $("#tree_table").treetable({ expandable: true });
    $("#tree_toggle").bind('click', function () {
      var target = $(this);
      if (target.attr("aria-pressed") == 'true') {
        target.html('<i class="fa fa-angle-down"></i> Expand All');
        $('#gdc-loading-icon').fadeIn(100);

        setTimeout(function () {
          $("#tree_table").find('a[title="Collapse"]').each(function () {
            $(this).trigger("click");
          });
          $('#gdc-loading-icon').fadeOut('fast');
        }, 1000);
      } else {
        target.html('<i class="fa fa-angle-up"></i>  Collapse All');
        $('#gdc-loading-icon').fadeIn(100);

        setTimeout(function () {
          $("#tree_table").find('a[title="Expand"]').each(function () {
            $(this).trigger("click");
          });
          $("#tree_table").find('a[title="Expand"]').each(function () {
            $(this).trigger("click");
          });
          $("#tree_table").find('a[title="Expand"]').each(function () {
            $(this).trigger("click");
          });
          $('#gdc-loading-icon').fadeOut('fast');
        }, 1000);
      }
    });
  }

  $('a.redirect').bind('click', function (event) {
    if (window.location.href.indexOf('https://docs.gdc.cancer.gov/') < 0) {
      event.preventDefault();
      var href = $(this).attr('href');
      window.open('https://docs.gdc.cancer.gov' + href, '_blank');
    }
  });

  $('#tab-values').bind('click', function () {
    var option = JSON.parse(localStorage.getItem('option'));
    option.activeTab = 0;
    localStorage.setItem('option', JSON.stringify(option));
  });

  $('#tab-properties').bind('click', function () {
    var option = JSON.parse(localStorage.getItem('option'));
    option.activeTab = 1;
    localStorage.setItem('option', JSON.stringify(option));
  });

  $('#tab-dictionary').bind('click', function () {
    var option = JSON.parse(localStorage.getItem('option'));
    option.activeTab = 2;
    localStorage.setItem('option', JSON.stringify(option));
  });

  var htmlShow = '';

  $('.show-more-less').click(function () {
    var target = $(this);

    var parentTable = $(this).parent().parent().parent();
    var targets = parentTable.find('.table__row--toggle');
    if (target.hasClass('more')) {
      target.removeClass('more');
      targets.slideToggle(350);
      target.html(htmlShow == '' ? '<i class="fa fa-angle-down"></i> Show More' : htmlShow);
    } else {
      htmlShow = target.html();
      target.addClass('more');
      targets.slideToggle(350).css({ display: 'flex' });
      target.html('<i class="fa fa-angle-up"></i> Show Less');
    }
  });

  $('.collapser').click(function () {
    var target = $(this);
    var parentTable = $(this).parent().parent().parent();

    var dataContainer = parentTable.find('#data-content');

    dataContainer.slideToggle(400, function () {
      if (dataContainer.is(":visible")) {
        target.html('<i class="fa fa-minus"></i>');
      } else {
        target.html('<i class="fa fa-plus"></i>');
      }
    });
  });

  $('.gdc-details').click(function () {
    var target = $(this);
    var parentTarget = $(this).parent();
    var gdcLinks = parentTarget.find('#gdc-links');
    gdcLinks.slideToggle(350);
  });

  var hiddenRows = $('#tree_table').find('.data-hide');
  $('#trs-checkbox').click(function () {
    if (this.checked) {
      hiddenRows.each(function () {
        $(this).removeClass('hide');
      });
    } else {
      hiddenRows.each(function () {
        $(this).addClass('hide');
      });
    }
  });

  $('.cde-suggest').click(function () {
    var alertSuggest = $('#alert-suggest');
    alertSuggest.removeClass('animated fadeInDownUp').css({ 'display': 'none' });
    var animationEnd = 'webkitAnimationEnd mozAnimationEnd MSAnimationEnd oanimationend animationend';
    alertSuggest.css({ 'display': 'block', 'top': __WEBPACK_IMPORTED_MODULE_4__shared__["a" /* default */].headerOffset() + 20 + 'px' }).addClass('animated fadeInDownUp').one(animationEnd, function () {
      alertSuggest.css({ 'display': 'none' });
    });
  });

  var windowEl = $(window);

  windowEl.resize(function () {
    var heightSlider = $('.navbar .container').height();
    var dialogs = $('#gdc_data, #gdc_syn_data, #compare_dialog, #caDSR_data, #compareGDC_dialog');
    dialogs.each(function () {
      var target = $(this).parent();
      if (target.offset().top < heightSlider) {
        target.css('top', heightSlider + 10 + "px");
      } else if (windowEl.width() < target.offset().left + target.width()) {
        target.css('left', windowEl.width() - target.width() - 10 + "px");
      }
    });
  });

  $('.tooltip-target').tooltip();
}

/***/ }),
/* 3 */
/***/ (function(module, exports, __webpack_require__) {

__webpack_require__(4);
module.exports = __webpack_require__(25);


/***/ }),
/* 4 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
Object.defineProperty(__webpack_exports__, "__esModule", { value: true });
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__search_bar___ = __webpack_require__(5);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__dialog___ = __webpack_require__(15);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__render__ = __webpack_require__(2);





$("#search").bind("click", __WEBPACK_IMPORTED_MODULE_0__search_bar___["a" /* default */].search);

$("#keywords").bind("keypress", __WEBPACK_IMPORTED_MODULE_0__search_bar___["a" /* default */].gotoSearch);

$("#keywords").bind("keydown", __WEBPACK_IMPORTED_MODULE_0__search_bar___["a" /* default */].selectSuggestion);

$("#keywords").bind("input", __WEBPACK_IMPORTED_MODULE_0__search_bar___["a" /* default */].suggest);

$(document).on('click', __WEBPACK_IMPORTED_MODULE_0__search_bar___["a" /* default */].removeBox);

var heightSlider = $('.navbar .container').height();
$('#docs-container').attr('style', 'margin-top: ' + (heightSlider - 54) + 'px !important');

$(window).resize(function () {
    heightSlider = $('.navbar .container').height();
    $('#docs-container').attr('style', 'margin-top: ' + (heightSlider - 54) + 'px !important');
});

function getGDCData(prop, target) {
    var uid = prop.replace(/@/g, '/');
    __WEBPACK_IMPORTED_MODULE_1__dialog___["a" /* default */].getGDCData(uid, target);
}

window.getGDCData = getGDCData;

function getGDCTerms(prop, targets) {
    var uid = prop.replace(/@/g, '/');
    __WEBPACK_IMPORTED_MODULE_1__dialog___["a" /* default */].getGDCTerms(uid, targets);
};

window.getGDCTerms = getGDCTerms;

function toCompare(prop) {
    var uid = prop.replace(/@/g, '/');
    __WEBPACK_IMPORTED_MODULE_1__dialog___["a" /* default */].toCompare(uid);
};

window.toCompare = toCompare;

function compare(gv) {
    if ($('#cp_input').val().trim() === '') {
        $('#cp_massage').css("display", "block");
        $("#cp_massage").removeClass();
        $('#cp_massage').addClass("compare-form__message");
        $('#cp_massage').html("Please type in user defined values.");
        return;
    } else {
        //compare and render
        $('#cp_massage').css("display", "none");
        $("#cp_massage").removeClass();
        $('#cp_massage').html("");
        $('#compare_form').css("display", "none");
        $('#compare_result').css("display", "block");

        var compare_dialog = $('#compare_dialog').parent().find('.ui-dialog-titlebar');

        var titleComponent = '<div class="checkbox ui-checkbox"><label class="checkbox__label checkbox__label--height"><input id="compare_filter" class="checkbox__input" type="checkbox" value=""><span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Case Sensitive</label>' + '<label class="checkbox__label checkbox__label--height"><input id="compare_unmatched" class="checkbox__input" type="checkbox" value=""><span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Hide Unmatched Values</label>';

        compare_dialog.append(titleComponent);

        var vs = $('#cp_input').val().split(/\n/);

        var opt = {};
        opt.sensitive = false;
        opt.unmatched = false;
        var table = generateCompareResult(vs, gv, opt);
        var html = '' //'<div id="cp_result_option"><div class="option-left"><input type="checkbox" id="compare_filter"> Case Sensitive</div><div class="option-right"><input type="checkbox" id="compare_unmatched"> Hide Unmatched Values</div></div>'
        + '<div id="cp_result_table" class="table__container table__container--margin-bottom">' + table + '</div>' + '<div id="cp_result_bottom" class="compare_result__bottom"><button id="back2Compare" class="btn btn-default compare_result__button">Back</button></div>';

        $('#compare_result').html(html);

        $('#compare_filter').bind('click', function () {
            var options = {};
            options.sensitive = $("#compare_filter").prop('checked');
            options.unmatched = $("#compare_unmatched").prop('checked');
            var table_new = generateCompareResult(vs, gv, options);
            $('#cp_result_table').html(table_new);
        });
        $('#compare_unmatched').bind('click', function () {
            var options = {};
            options.sensitive = $("#compare_filter").prop('checked');
            options.unmatched = $("#compare_unmatched").prop('checked');
            var table_new = generateCompareResult(vs, gv, options);
            $('#cp_result_table').html(table_new);
        });
        $('#back2Compare').bind('click', function () {
            $('#compare_result').html("");
            $('#compare_result').css("display", "none");
            $('#compare_form').css("display", "block");
            compare_dialog.find('.ui-checkbox').remove();
        });
    }
};

window.compare = compare;

function generateCompareResult(fromV, toV, option) {
    var v_lowercase = [],
        v_matched = [];
    if (option.sensitive) {
        toV.forEach(function (v) {
            v_lowercase.push(v.trim());
        });
    } else {
        toV.forEach(function (v) {
            v_lowercase.push(v.trim().toLowerCase());
        });
    }

    var table = '<div class="table__thead row">' + '<div class="table__th col-xs-6">User Defined Values</div>' + '<div class="table__th col-xs-6">Matched GDC Values</div>' + '</div>' + '<div class="table__body row">' + '<div class="col-xs-12">';

    fromV.forEach(function (v) {
        var tmp = $.trim(v);
        if (tmp === '') {
            return;
        }
        var text = '';
        var idx = option.sensitive ? v_lowercase.indexOf(tmp) : v_lowercase.indexOf(tmp.toLowerCase());
        if (idx >= 0) {
            text = toV[idx];
            v_matched.push(idx);
        }
        if (text === '') {
            text = '<div style="color:red;">--</div>';
            //table += '<tr class="data-table-row"><td align="left">'+v+'</td><td align="left">'+text+'</td></tr>';
            table += '<div class="table__row row">' + '<div class="table__td table__td--slim col-xs-6">' + v + '</div>' + '<div class="table__td table__td--slim col-xs-6">' + text + '</div>' + '</div>';
        } else {
            //table += '<tr class="data-table-row"><td align="left">'+v+'</td><td align="left"><b>'+(idx+1)+'.</b>'+text+'</td></tr>';
            table += '<div class="table__row row">' + '<div class="table__td table__td--slim col-xs-6">' + v + '</div>' + '<div class="table__td table__td--slim col-xs-6">' + text + '</div>' + '</div>';
        }
    });
    for (var i = 0; i < toV.length; i++) {
        if (v_matched.indexOf(i) >= 0) {
            continue;
        }
        //table += '<tr class="data-table-row '+(option.unmatched ? 'row-undisplay' : '')+'"><td align="left"><div style="color:red;">--</div></td><td align="left"><b>'+(i+1)+'.</b>'+toV[i]+'</td></tr>';
        table += '<div class="table__row row ' + (option.unmatched ? 'table__row--undisplay' : '') + '">' + '<div class="table__td table__td--slim col-xs-6"><div style="color:red;">--</div></div>' + '<div class="table__td table__td--slim col-xs-6">' + toV[i] + '</div>' + '</div>';
    }
    table += '</div></div>';
    //table += "</tbody></table>";
    return table;
};

window.generateCompareResult = generateCompareResult;

function generateCompareGDCResult(fromV, toV, option) {
    var v_lowercase = [],
        v_matched = [];
    var from_num = 0;
    if (option.sensitive) {
        toV.forEach(function (v) {
            v_lowercase.push(v.trim());
        });
    } else {
        toV.forEach(function (v) {
            v_lowercase.push(v.trim().toLowerCase());
        });
    }

    var table = '<div class="table__thead row">' + '<div class="table__th col-xs-6">GDC Values</div>' + '<div class="table__th col-xs-6">Matched caDSR Values</div>' + '</div>' + '<div class="table__body row">' + '<div class="col-xs-12">';

    fromV.forEach(function (v) {
        var tmp = $.trim(v);
        if (tmp === '') {
            return;
        }
        var text = '';
        var idx = option.sensitive ? v_lowercase.indexOf(tmp) : v_lowercase.indexOf(tmp.toLowerCase());
        if (idx >= 0) {
            text = toV[idx];
            v_matched.push(idx);
        }
        if (text === '') {
            text = '<div style="color:red;">--</div>';
            table += '<div class="table__row row">' + '<div class="table__td table__td--slim col-xs-6">' + v + '</div>' + '<div class="table__td table__td--slim col-xs-6">' + text + '</div>' + '</div>';
        } else {
            table += '<div class="table__row row">' + '<div class="table__td table__td--slim col-xs-6">' + v + '</div>' + '<div class="table__td table__td--slim col-xs-6">' + text + '</div>' + '</div>';
        }
    });
    for (var i = 0; i < toV.length; i++) {
        if (v_matched.indexOf(i) >= 0) {
            continue;
        }
        table += '<div class="table__row row ' + (option.unmatched ? 'table__row--undisplay' : '') + '">' + '<div class="table__td  table__td--slim col-xs-6"><div style="color:red;">--</div></div>' + '<div class="table__td table__td--slim col-xs-6">' + toV[i] + '</div>' + '</div>';
    }
    table += '</div></div>';

    return table;
};

window.generateCompareGDCResult = generateCompareGDCResult;

function getCDEData(cdeId, targets) {
    __WEBPACK_IMPORTED_MODULE_1__dialog___["a" /* default */].getCDEData(cdeId, targets);
}

window.getCDEData = getCDEData;

function compareGDC(prop, cdeId) {
    var uid = prop.replace(/@/g, '/');
    __WEBPACK_IMPORTED_MODULE_1__dialog___["a" /* default */].compareGDC(uid, cdeId);
};

window.compareGDC = compareGDC;

function getNCITDetails(uid) {
    __WEBPACK_IMPORTED_MODULE_1__dialog___["a" /* default */].getNCITDetails(uid);
}

window.getNCITDetails = getNCITDetails;

//find the word with the first character capitalized
function findWord(words) {
    var word = "";
    if (words.length == 1) {
        return words[0];
    }
    words.forEach(function (w) {
        if (word !== "") {
            return;
        }
        var idx_space = w.indexOf(" ");
        var idx_comma = w.indexOf(",");
        if (idx_space == -1 && idx_comma == -1) {
            if (/^[A-Z][a-z0-9]{0,}$/.test(w)) {
                word = w;
            }
        } else if (idx_space !== -1 && idx_comma == -1) {
            if (/^[A-Z][a-z0-9]{0,}$/.test(w.substr(0, idx_space))) {
                word = w;
            }
        } else if (idx_space == -1 && idx_comma !== -1) {
            if (/^[A-Z][a-z0-9]{0,}$/.test(w.substr(0, idx_comma))) {
                word = w;
            }
        } else {
            if (idx_comma > idx_space) {
                if (/^[A-Z][a-z0-9]{0,}$/.test(w.substr(0, idx_space))) {
                    word = w;
                }
            } else {
                if (/^[A-Z][a-z0-9]{0,}$/.test(w.substr(0, idx_comma))) {
                    word = w;
                }
            }
        }
    });
    if (word == "") {
        word = words[0];
    }
    return word;
};

window.findWord = findWord;

$(function () {

    $('#body a[href^="http"]').each(function () {
        var anchor = $(this);
        anchor.removeClass('external-link');
        anchor.html($.trim(anchor[0].innerText));
    });

    if (localStorage.hasOwnProperty('keyword') || localStorage.hasOwnProperty('option') || localStorage.hasOwnProperty('items')) {

        $('#gdc-loading-icon').fadeIn(100);

        setTimeout(function () {

            var keyword = localStorage.getItem('keyword');
            var option = JSON.parse(localStorage.getItem('option'));
            var items = JSON.parse(localStorage.getItem('items'));

            if (keyword != null || option != null || items != null) {

                $("#keywords").val(keyword);

                if (option.match != 'partial') {
                    $("#i_syn").prop('checked', true);
                }
                if (option.desc != false) {
                    $("#i_desc").prop('checked', true);
                }
                if (option.syn != false) {
                    $("#i_syn").prop('checked', true);
                }

                Object(__WEBPACK_IMPORTED_MODULE_2__render__["a" /* default */])(keyword, option, items);

                $('#gdc-loading-icon').fadeOut('fast');
            }
        }, 100);
    }
});

/***/ }),
/* 5 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__api__ = __webpack_require__(0);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__render__ = __webpack_require__(2);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__view__ = __webpack_require__(14);




var displayBoxIndex = -1;
var activeTab = 0;

var func = {
    search: function search() {
        var keyword = $("#keywords").val();
        var option = {};

        if (keyword == "") {
            option.error = true;
            $('#keywords').addClass('search-bar__input--has-error');
            Object(__WEBPACK_IMPORTED_MODULE_1__render__["a" /* default */])(keyword, option, []);
            return;
        }

        keyword = keyword.toLowerCase();
        //get selected tab
        var count = 0;
        $("li[role='presentation']").each(function () {
            if ($(this).hasClass("active")) {
                activeTab = count;
            }
            count++;
        });

        option.desc = $("#i_desc").prop('checked');
        option.syn = $("#i_syn").prop('checked');
        option.match = $("#i_ematch").prop('checked') ? "exact" : "partial";
        option.activeTab = option.desc ? 1 : activeTab;
        $("#suggestBox").css("display", "none");
        displayBoxIndex = -1;
        //todo:show progress bar
        $('#gdc-loading-icon').fadeIn(100);
        __WEBPACK_IMPORTED_MODULE_0__api__["a" /* default */].searchAll(keyword, option, function (keyword, option, items) {

            //Save the data in localStorage
            localStorage.setItem('keyword', keyword);
            localStorage.setItem('option', JSON.stringify(option));
            localStorage.setItem('items', JSON.stringify(items));

            Object(__WEBPACK_IMPORTED_MODULE_1__render__["a" /* default */])(keyword, option, items);
            //todo: close progress bar
            $('#gdc-loading-icon').fadeOut('fast');
        });
    },
    gotoSearch: function gotoSearch(e) {
        if (e.keyCode == 13) {
            e.preventDefault();
        }
        if (e.keyCode == 13 && $("#suggestBox .selected").length !== 0) {
            var t = $("#suggestBox .selected").text();
            $("#keywords").val(t.substr(0, t.length - 1));
            $("#search").trigger("click");
        } else if (e.keyCode == 13) {
            $("#search").trigger("click");
        }
    },
    selectSuggestion: function selectSuggestion(e) {
        if ((e.keyCode == 40 || e.keyCode == 38) && $(this).val().trim() !== "" && document.getElementById("suggestBox").style.display !== "none") {
            e.preventDefault();
            //focus to the first element

            displayBoxIndex += e.keyCode == 40 ? 1 : -1;
            var oBoxCollection = $("#suggestBox").find("div");
            if (displayBoxIndex >= oBoxCollection.length) displayBoxIndex = 0;
            if (displayBoxIndex < 0) displayBoxIndex = oBoxCollection.length - 1;
            var cssClass = "selected";
            oBoxCollection.removeClass(cssClass).eq(displayBoxIndex).addClass(cssClass);
        }
    },
    suggest: function suggest() {
        var area = document.getElementById("suggestBox");

        if ($("#keywords").hasClass('search-bar__input--has-error')) {
            $("#keywords").removeClass('search-bar__input--has-error');
        }

        if ($(this).val().trim() === '') {
            area.style.display = "none";
            displayBoxIndex = -1;
            area.innerHTML = "";
            return;
        }

        __WEBPACK_IMPORTED_MODULE_0__api__["a" /* default */].suggest($(this).val(), function (result) {
            if (result.length === 0) {
                area.style.display = "none";
                displayBoxIndex = -1;
                area.innerHTML = "";
                return;
            }

            area.style.display = "block";
            var html = $.templates(__WEBPACK_IMPORTED_MODULE_2__view__["a" /* default */]).render({ options: result });;
            displayBoxIndex = -1;
            area.innerHTML = html;
            area.onclick = function (e) {
                var t = $(e.target).text();
                $("#keywords").val(t);
                $("#keywords").focus();
            };
        });
    },
    removeBox: function removeBox(e) {
        if ($(e.target) != $("#suggestBox")) {
            $("#suggestBox").css("display", "none");
            displayBoxIndex = -1;
        }
    }
};

/* harmony default export */ __webpack_exports__["a"] = (func);

/***/ }),
/* 6 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(7);


var func = {
    render: function render(items) {
        //data preprocessing
        //current category
        var c_c = "";
        //current node
        var c_n = "";
        //prefix for property and value id
        var count = 0;
        //data generated
        var trs = [];

        //category row
        var c = {};
        //node row
        var n = {};
        //final result
        var result = {};
        result.len = 0;

        items.forEach(function (item) {
            var hl = item.highlight;
            var source = item._source;
            var enum_s = "enum.s" in hl || "enum.s.have" in hl ? hl['enum.s'] || hl["enum.s.have"] : [];
            var enum_n = "enum.n" in hl || "enum.n.have" in hl ? hl["enum.n"] || hl["enum.n.have"] : [];
            var cde_n = "cde_pv.n" in hl || "cde_pv.n.have" in hl ? hl["cde_pv.n"] || hl["cde_pv.n.have"] : [];
            var cde_s = "cde_pv.ss.s" in hl || "cde_pv.ss.s.have" in hl ? hl["cde_pv.ss.s"] || hl["cde_pv.ss.s.have"] : [];
            var arr_enum_s = [];
            var arr_enum_n = [];
            var arr_cde_n = [];
            var arr_cde_s = [];
            var matched_pv = [];
            enum_s.forEach(function (s) {
                var tmp = s.replace(/<b>/g, "").replace(/<\/b>/g, "");
                arr_enum_s.push(tmp);
            });
            enum_n.forEach(function (n) {
                var tmp = n.replace(/<b>/g, "").replace(/<\/b>/g, "");
                arr_enum_n.push(tmp);
            });
            cde_n.forEach(function (pn) {
                var tmp = pn.replace(/<b>/g, "").replace(/<\/b>/g, "");
                arr_cde_n.push(tmp);
            });
            cde_s.forEach(function (ps) {
                var tmp = ps.replace(/<b>/g, "").replace(/<\/b>/g, "");
                arr_cde_s.push(tmp);
            });

            if (source.cde_pv !== undefined && source.cde_pv.length > 0) {
                source.cde_pv.forEach(function (pv) {
                    var exist = false;
                    if (pv.ss !== undefined && pv.ss.length > 0) {
                        pv.ss.forEach(function (ss) {
                            ss.s.forEach(function (s) {
                                if (arr_cde_s.indexOf(s) !== -1) {
                                    exist = true;
                                }
                            });
                        });
                    }
                    exist = exist || arr_cde_n.indexOf(pv.n) >= 0;
                    if (exist) {
                        matched_pv.push(pv.n.toLowerCase());
                    }
                });
            }

            if (source.category != c_c) {
                //put category to tree table
                count++;
                c_c = source.category;
                c = {};
                c.id = c_c;
                c.title = c_c;
                c.desc = "";
                c.data_tt_id = count + "_" + c.id;
                c.data_tt_parent_id = "--";
                c.type = "category";
                c.node = "branch";
                c.exist = true;
                c.len = 0;
                trs.push(c);
            }
            if (source.node != c_n) {
                //put node to tree table
                count++;
                c_n = source.node;
                n = {};
                //link id
                n.l_id = source.node;
                n.id = source.node;
                n.title = source.n_title;
                n.desc = source.n_desc;
                n.data_tt_id = count + "_" + n.id;
                n.data_tt_parent_id = c.data_tt_id;
                n.type = "folder";
                n.node = "branch";
                n.exist = true;
                n.len = 0;
                trs.push(n);
            }
            //put property to tree table
            var p = {};
            //calculate if property itself got matched;
            var count_p = 0;
            //calculate how many values got matched in this property;
            var count_v = 0;
            //calculate how many synonyms got matched in this property
            var count_s = 0;
            count++;
            p.id = count + "_" + source.name;
            //link id
            p.l_id = source.name;
            p.parent_l_id = n.l_id;
            //may have highlighted terms in p.title and p.desc
            p.title = "name" in hl || "name.have" in hl ? hl["name"] || hl["name.have"] : source.name;
            p.desc = "desc" in hl ? hl["desc"] : source.desc;
            p.data_tt_id = p.id;
            p.data_tt_parent_id = n.data_tt_id;
            p.type = "property";
            p.exist = true;
            if ("name" in hl || "name.have" in hl || "desc" in hl) {
                count_p = 1;
            }
            //put value to tree table
            if (source.enum != undefined) {
                if (enum_n.length == 0 && enum_s.length == 0 && matched_pv.length == 0) {
                    //if no values show in the values tab
                    p.node = "";
                    trs.push(p);
                } else {
                    p.node = "branch";
                    trs.push(p);
                    //show values, need to highlight if necessary
                    var list = [];
                    if ("enum.n" in hl || "enum.n.have" in hl) {
                        list = hl["enum.n"] || hl["enum.n.have"];
                    }
                    var enums = {};
                    list.forEach(function (em) {
                        var e = em.replace(/<b>/g, "").replace(/<\/b>/g, "");
                        enums[e] = em;
                    });
                    var values = source.enum;
                    var tmp_trs = [];
                    values.forEach(function (v) {
                        count++;
                        var e = {};
                        e.id = count + "_" + v.n;
                        e.exist = false;

                        var idx = matched_pv.indexOf(v.n.toLowerCase());
                        if (idx !== -1) {
                            count_s--;
                            e.exist = true;
                        } else {
                            if (arr_enum_n.indexOf(v.n) !== -1) {
                                e.exist = true;
                            }

                            if (v.s !== undefined && e.exist != true) {
                                v.s.forEach(function (syn) {
                                    if (arr_enum_s.indexOf(syn) !== -1) {
                                        e.exist = true;
                                    }
                                });
                            }
                        }

                        if (e.exist) {
                            count_v++;
                        }
                        //may be highlighted
                        e.title = v.n in enums ? enums[v.n] : v.n;
                        e.desc = "";
                        e.data_tt_id = e.id;
                        e.data_tt_parent_id = p.id;
                        e.type = "value";
                        e.node = "leaf";
                        tmp_trs.push(e);
                    });
                    if (count_v == 0) {
                        p.node = "";
                    } else {
                        tmp_trs.forEach(function (tt) {
                            trs.push(tt);
                        });
                    }

                    count_s += matched_pv.length;
                }
            } else {
                if (matched_pv.length > 0) {
                    //matched on cde
                    p.node = "branch";
                    trs.push(p);
                    //show caDSR reference
                    count++;
                    var l = {};
                    l.id = count + "_l";
                    l.l_id = source.cde.id;
                    l.l_type = "cde";
                    l.url = source.cde.url;
                    l.desc = "";
                    l.data_tt_id = l.id;
                    l.data_tt_parent_id = p.id;
                    l.type = "link";
                    l.node = "leaf";
                    trs.push(l);

                    count_s = matched_pv.length;
                } else {
                    //matched on property name or description
                    p.node = "";
                    trs.push(p);
                }
            }

            //save and calculate the count of matched element in this property
            p.len = count_v + count_s;
            c.len += p.len + count_p;
            n.len += p.len + count_p;
            result.len += p.len + count_p;
        });

        var offset = $('#root').offset().top;
        var h = window.innerHeight - offset - 305;
        h = h < 430 ? 430 : h;
        var html = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */]).render({ mh: h, trs: trs });

        result.html = html;
        return result;
    }
};

/* harmony default export */ __webpack_exports__["a"] = (func);

/***/ }),
/* 7 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

var tmpl = '<div class="container table__container">' + '<div class="table__thead row">' + '<div class="col-xs-4"><div class="table__th">Name</div></div>' + '<div class="col-xs-4"><div class="table__th">Description</div></div>' + '<div class="col-xs-4"><div class="table__th table__th--right">' + '<div class="checkbox checkbox-th"><label class="checkbox__label">' + '<input class="checkbox__input" id="trs-checkbox" type="checkbox" value="">' + '<span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Show all values' + '</label></div>' + '<button type="button" id="tree_toggle" class="btn btn-th" data-toggle="button" aria-pressed="false" autocomplete="off">' + '<i class="btn-th__icon fa fa-angle-down"></i> Expand All' + '</button>' + '</div></div>' + '</div>' + '<div class="table__body table__body--overflow row" style="max-height: {{:mh}}px;">' + '<table class="treetable table" id="tree_table">' + '<tbody>' + '{{for trs}}' + '<tr key="{{:id}}" data-tt-id="{{:data_tt_id}}" data-tt-parent-id="{{:data_tt_parent_id}}" class="data-table-row {{:node}} {{if exist != true && type == "value"}}data-hide hide{{/if}}">' + '<td width="33%">' + '<span class="{{:type}} title table__td--word-break" style="display:inline-block; width: 70%;">' + '{{if type == "category"}}' + '<a class="redirect" href="/Data_Dictionary/viewer/#?view=table-entity-list&anchor={{:id}}">{{:title}} {{if len !== undefined}}({{:len}}){{/if}}</a>' + '{{else type == "folder"}}' + '<a class="redirect" href="/Data_Dictionary/viewer/#?view=table-definition-view&id={{:l_id}}">{{:title}} {{if len !== undefined}}({{:len}}){{/if}}</a>' + '{{else type == "property"}}' + '<a class="redirect" href="/Data_Dictionary/viewer/#?view=table-definition-view&id={{:parent_l_id}}&anchor={{:l_id}}">{{:title}} {{if len !== undefined}}({{:len}}){{/if}}</a>' + '{{else type == "link"}}' + '{{if l_type == "cde"}}' + 'No values in GDC, reference values in <a href="javascript:getCDEData(\'{{:l_id}}\');" class="table-td-link">caDSR</a>' + '{{else}}' + 'No values in GDC, concept referenced in <a target="_blank" href="{{:url}}" class="table-td-link">NCIt</a>' + '{{/if}}' + '{{else}}' + '{{:title}} {{if len && len !== 0}}({{:len}}){{/if}}' + '{{/if}}' + '</td>' + '<td width="66%">' + '{{:desc}}' + '</td>' + '</tr>' + '{{/for}}' + '</tbody>' + '</table>' + '</div>' + '</div>';

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 8 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(9);


var func = {
  render: function render(items) {
    //data preprocessing
    var props = [];
    items.forEach(function (item) {
      var hl = item.highlight;
      var source = item._source;
      if ("name" in hl || "name.have" in hl || "desc" in hl) {
        var prop = {};
        prop.nm = "name" in hl || "name.have" in hl ? hl["name"] || hl["name.have"] : source.name;
        prop.nd = source.node;
        prop.ct = source.category;
        prop.desc = "desc" in hl ? hl["desc"] : source.desc;
        prop.local = source.enum == undefined ? false : true;
        prop.syn = false;
        if (source.enum !== undefined) {
          //check if synonyms exists
          source.enum.forEach(function (em) {
            if (prop.syn) return;

            if (em.n_c !== undefined) {
              prop.syn = true;
            }
          });
        }
        prop.ref = source.name + "@" + source.node + "@" + source.category;
        prop.cdeId = source.cde !== undefined ? source.cde.id : "";
        prop.cdeUrl = source.cde !== undefined ? source.cde.url : "";
        prop.cdeLen = source.cde_pv == undefined || source.cde_pv.length == 0 ? false : true;
        prop.type = Array.isArray(source.type) ? source.type[0] : source.type;
        if (source.cde !== undefined && source.cde.dt !== undefined) {
          prop.type = source.cde.dt;
        }
        if (prop.type) {
          prop.type = prop.type.toLowerCase();
        }

        props.push(prop);
      }
    });
    var html = "";
    if (props.length == 0) {
      var keyword = $("#keywords").val();
      html = '<div class="indicator">Sorry, no results found for kerword: <span class="indicator__term">' + keyword + '</span></div>';
    } else {
      var offset = $('#root').offset().top;
      var h = window.innerHeight - offset - 300;
      h = h < 430 ? 430 : h;
      html = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */]).render({ mh: h, props: props });
    }

    var result = {};
    result.len = props.length;
    result.html = html;
    return result;
  }
};

/* harmony default export */ __webpack_exports__["a"] = (func);

/***/ }),
/* 9 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

var tmpl = '<div class="container table__container"><div class="table__thead row table__thead--padding-right">' + '<div class="table__th col-xs-2">Category / Node</div>' + '<div class="table__th col-xs-2">Property</div>' + '<div class="table__th col-xs-4">Description</div>' + '<div class="table__th col-xs-2">GDC Property Values</div>' + '<div class="table__th col-xs-2">caDSR CDE Reference</div>' + '</div>' + '<div class="table__body table__body--overflow row" style="max-height: {{:mh}}px;"><div class="col-xs-12">' + '{{for props}}' + '<div class="table__row row">' + '<div class="table__td col-xs-2">' + '<a class="redirect" href="/Data_Dictionary/viewer/#?view=table-entity-list&anchor={{:ct}}">{{:ct}}</a>' + '<ul><li class="table__td--word-break">' + '<a class="redirect" href="/Data_Dictionary/viewer/#?view=table-definition-view&id={{:nd}}">{{:nd}}</a>' + '</li></ul>' + '</div>' + '<div class="table__td col-xs-2 table__td--word-break">' + '<a class="redirect" href="/Data_Dictionary/viewer/#?view=table-definition-view&id={{:nd}}&anchor={{:nm}}">{{:nm}}</a>' + '</div>' + '<div class="table__td col-xs-4">{{:desc}}</div>' + '<div class="table__td col-xs-2">' + '{{if local}}' + '<a href="javascript:getGDCData(\'{{:ref}}\',null);">See All Values</a>' + '<br><a href="javascript:toCompare(\'{{:ref}}\');"> Compare with User List</a>' + '{{if syn}}' + '<br><a href="javascript:getGDCTerms(\'{{:ref}}\',null);">See All Terms</a>' + '{{else}}' + '{{/if}}' + '{{else}}' + 'type: {{:type}}' + '{{/if}}' + '</div>' + '<div class="table__td col-xs-2">' + '{{if cdeId == ""}}' + '' + '{{else}}' + '<a class="table-td-link" href="{{:cdeUrl}}" target="_blank">CDE</a>' + '{{if local && cdeLen}}' + ', <a class="table-td-link" href="javascript:getCDEData(\'{{:cdeId}}\',null);">Values</a>, <a class="table-td-link" href="javascript:compareGDC(\'{{:ref}}\',\'{{:cdeId}}\');"> Compare with GDC</a>' + '{{else cdeLen}}' + ', <a class="table-td-link" href="javascript:getCDEData(\'{{:cdeId}}\',null);">Values</a>' + '{{else}}' + '' + '{{/if}}' + '{{/if}}' + '</div>' + '</div>' + '{{/for}}</div></div></div>';

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 10 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(11);


var func = {
  render: function render(items) {
    //data preprocessing

    var values = [];
    var len = 0;
    items.forEach(function (item) {
      var hl = item.highlight;
      if (hl["enum.n"] == undefined && hl["enum.n.have"] == undefined && hl["enum.s"] == undefined && hl["enum.s.have"] == undefined && hl["cde_pv.n"] == undefined && hl["cde_pv.n.have"] == undefined && hl["cde_pv.ss.s"] == undefined && hl["cde_pv.ss.s.have"] == undefined && hl["enum.i_c.c"] == undefined && hl["enum.i_c.have"] == undefined) {
        return;
      }
      var source = item._source;
      var dict_enum_n = {};
      var dict_enum_s = {};
      var dict_cde_n = {};
      var dict_cde_s = {};
      var arr_enum_c = [];
      var arr_enum_c_have = [];
      //each row in the values tab will be put into values
      var row = {};
      row.category = source.category;
      row.node = source.node;
      row.name = source.name;
      row.local = source.enum == undefined ? false : true;
      row.syn = false;
      if (source.enum !== undefined) {
        //check if synonyms exists
        source.enum.forEach(function (em) {
          if (row.syn) return;

          if (em.n_c !== undefined) {
            row.syn = true;
          }
        });
      }
      row.ref = source.name + "@" + source.node + "@" + source.category;
      row.cdeId = source.cde !== undefined ? source.cde.id : "";
      row.cdeUrl = source.cde !== undefined ? source.cde.url : "";
      row.cdeLen = source.cde_pv == undefined || source.cde_pv.length == 0 ? false : true;
      //value informations in the subtable
      row.vs = [];
      row.tgts_enum_n = ""; //added
      row.tgts_cde_n = "";
      var enum_n = "enum.n" in hl || "enum.n.have" in hl ? hl["enum.n"] || hl["enum.n.have"] : [];
      var enum_s = "enum.s" in hl || "enum.s.have" in hl ? hl['enum.s'] || hl["enum.s.have"] : [];
      var cde_n = "cde_pv.n" in hl || "cde_pv.n.have" in hl ? hl["cde_pv.n"] || hl["cde_pv.n.have"] : [];
      var cde_s = "cde_pv.ss.s" in hl || "cde_pv.ss.s.have" in hl ? hl["cde_pv.ss.s"] || hl["cde_pv.ss.s.have"] : [];
      var enum_c = "enum.i_c.c" in hl ? hl["enum.i_c.c"] : [];
      var enum_c_have = "enum.i_c.have" in hl ? hl["enum.i_c.have"] : [];
      enum_n.forEach(function (n) {
        var tmp = n.replace(/<b>/g, "").replace(/<\/b>/g, "");
        dict_enum_n[tmp] = n;
      });
      enum_s.forEach(function (s) {
        var tmp = s.replace(/<b>/g, "").replace(/<\/b>/g, "");
        dict_enum_s[tmp] = s;
      });
      cde_n.forEach(function (pn) {
        var tmp = pn.replace(/<b>/g, "").replace(/<\/b>/g, "");
        dict_cde_n[tmp] = pn;
      });
      cde_s.forEach(function (ps) {
        var tmp = ps.replace(/<b>/g, "").replace(/<\/b>/g, "");
        dict_cde_s[tmp] = ps;
      });
      enum_c.forEach(function (c) {
        var tmp = c.replace(/<b>/g, "").replace(/<\/b>/g, "");
        if (arr_enum_c.indexOf(tmp) == -1) {
          arr_enum_c.push(tmp);
        }
      });
      enum_c_have.forEach(function (ch) {
        var tmp = ch.replace(/<b>/g, "").replace(/<\/b>/g, "");
        if (arr_enum_c_have.indexOf(tmp) == -1) {
          arr_enum_c_have.push(tmp);
        }
      });

      //check if there are any matches in the cde synonyms
      var matched_pv = {};
      if (source.cde_pv !== undefined && source.cde_pv.length > 0) {
        source.cde_pv.forEach(function (pv) {
          var exist = false;
          var tmp_ss = [];
          if (pv.ss !== undefined && pv.ss.length > 0) {
            pv.ss.forEach(function (ss) {
              var tmp_s = [];
              var tmp_s_h = [];
              //remove duplicate
              var cache = {};
              ss.s.forEach(function (s) {
                var lc = s.trim().toLowerCase();
                if (!(lc in cache)) {
                  cache[lc] = [];
                }
                cache[lc].push(s);
              });
              for (var idx in cache) {
                //find the term with the first character capitalized
                var word = findWord(cache[idx]);
                tmp_s.push(word);
              }
              tmp_s.forEach(function (s) {
                if (s in dict_cde_s) {
                  exist = true;
                  tmp_s_h.push(dict_cde_s[s]);
                } else {
                  tmp_s_h.push(s);
                }
              });
              tmp_ss.push({
                c: ss.c,
                s: tmp_s_h
              });
            });
          }
          exist = exist || pv.n in dict_cde_n;
          if (exist) {
            //matched_pv[pv.n.toLowerCase()] = tmp_ss;
            matched_pv[pv.n.toLowerCase()] = {
              "pv": pv.n in dict_cde_n ? dict_cde_n[pv.n] : pv.n,
              "pvm": pv.m,
              "ss": tmp_ss
            };
            pv.n = pv.n.replace(/\'/g, '^');
            row.tgts_cde_n += pv.n + "#";
          }
        });
      }

      if (source.enum) {
        source.enum.forEach(function (em) {
          //check if there are any matches in local synonyms
          var exist = false;
          var tmp_s = [];
          var t_s = [];
          if (em.s) {
            //remove depulicates in local synonyms
            var cache = {};
            em.s.forEach(function (s) {
              var lc = s.trim().toLowerCase();
              if (!(lc in cache)) {
                cache[lc] = [];
              }
              cache[lc].push(s);
            });
            for (var idx in cache) {
              //find the term with the first character capitalized
              var word = findWord(cache[idx]);
              t_s.push(word);
            }
            t_s.forEach(function (s) {
              if (s in dict_enum_s) {
                exist = true;
                tmp_s.push(dict_enum_s[s]);
              } else {
                tmp_s.push(s);
              }
            });
          }
          //value to be put into the subtable
          var v = {};
          if (exist) {
            //check if there is a match to the value name
            if (em.n in dict_enum_n) {
              v.n = dict_enum_n[em.n];
            } else {
              v.n = em.n;
            }
            v.ref = row.ref;
            v.n_c = em.n_c;
            v.s = tmp_s;
          } else {
            if (em.n in dict_enum_n) {
              v.n = dict_enum_n[em.n];
              v.ref = row.ref;
              v.n_c = em.n_c;
              //v.s = em.s;
              v.s = tmp_s;
            }
          }

          //check if it contains icd-0-3 codes.
          if (em.i_c !== undefined) {
            if (arr_enum_c.indexOf(em.i_c.c) >= 0) {
              v.i_c = "<b>" + em.i_c.c + "</b>";
              if (v.n == undefined) {
                v.n = em.n;
                v.ref = row.ref;
                v.n_c = em.n_c;
                //v.s = em.s;
                v.s = tmp_s;
              }
            } else {
              var has = false;
              em.i_c.have.forEach(function (ch) {
                if (has) return;
                if (arr_enum_c_have.indexOf(ch) >= 0) {
                  has = true;
                }
              });
              if (has) {
                v.i_c = "<b>" + em.i_c.c + "</b>";
                if (v.n == undefined) {
                  v.n = em.n;
                  v.ref = row.ref;
                  v.n_c = em.n_c;
                  //v.s = em.s;
                  v.s = tmp_s;
                }
              } else {
                v.i_c = em.i_c.c;
              }
            }
          }

          var lc = em.n.toLowerCase();
          if (lc in matched_pv) {
            if (v.n == undefined) {
              v.n = em.n;
              v.ref = row.ref;
              v.n_c = em.n_c;
              //v.s = em.s;
              v.s = tmp_s;
            }

            v.cde_s = matched_pv[lc].ss;
            if (v.cde_s.length) {
              v.cde_pv = matched_pv[lc].pv;
              v.cde_pvm = matched_pv[lc].pvm;
            }
            delete matched_pv[lc];
          } else {
            v.cde_s = [];
          }

          if (v.n !== undefined) {
            var tmp = v.n.replace(/<b>/g, "").replace(/<\/b>/g, "");
            row.tgts_enum_n += tmp + "#";
            row.vs.push(v);
          }
        });

        //add the rest of the matched cde_pvs to the subtables
        for (var idx in matched_pv) {
          var v = {};
          v.n = "no match";
          v.ref = row.ref;
          v.n_c = "";
          v.s = [];
          v.cde_s = matched_pv[idx].ss;
          if (v.cde_s.length) {
            v.cde_pv = matched_pv[idx].pv;
            v.cde_pvm = matched_pv[idx].pvm;
          }
          row.vs.push(v);
        }
        len += row.vs.length;

        //reformat the icd-o-3 code data
        if (row.vs) {
          var temp_i_c = [];
          var new_vs = [];
          row.vs.forEach(function (item) {
            if (item.i_c === undefined) {
              return;
            }
            temp_i_c.push(item.i_c.replace(/<b>/g, "").replace(/<\/b>/g, ""));
          });
          var results = [];
          temp_i_c.forEach(function (d) {
            if (results.indexOf(d) == -1) {
              results.push(d);
            }
          });

          if (results) {
            results.forEach(function (item) {

              var tmp_data = {
                i_c: {},
                n: [],
                ref: {},
                n_t: [],
                temp_n_c: []
              };

              row.vs.forEach(function (value) {
                if (value.i_c.replace(/<b>/g, "").replace(/<\/b>/g, "") == item) {
                  var temp_nt = {
                    n_c: {},
                    s: []
                  };
                  tmp_data.i_c = value.i_c;
                  tmp_data.cde_s = value.cde_s;
                  tmp_data.ref = value.ref;
                  tmp_data.n.push(value.n);
                  //tmp_data.ref.push(value.ref);
                  if (value.n_c && tmp_data.temp_n_c.indexOf(value.n_c) == -1) {
                    tmp_data.temp_n_c.push(value.n_c);
                    temp_nt.n_c = value.n_c;
                    value.s.forEach(function (syn) {
                      temp_nt.s.push(syn);
                    });

                    tmp_data.n_t.push(temp_nt);
                  }
                }
              });
              new_vs.push(tmp_data);
            });
          }
          if (new_vs.length !== 0) {
            row.vs = new_vs;
          }
        }
      }

      values.push(row);
    });
    var html = "";
    if (values.length == 0) {
      var keyword = $("#keywords").val();
      html = '<div class="indicator">Sorry, no results found for kerword: <span class="indicator__term">' + keyword + '</span></div>';
    } else {
      var offset = $('#root').offset().top;
      var h = window.innerHeight - offset - 300;
      h = h < 430 ? 430 : h;
      html = $.templates({
        markup: __WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */],
        allowCode: true
      }).render({
        mh: h,
        values: values
      });
    }
    var result = {};
    result.len = len;
    result.html = html;
    return result;
  }
};

/* harmony default export */ __webpack_exports__["a"] = (func);

/***/ }),
/* 11 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

var tmpl = '<div class="container table__container"><div class="table__thead row">' + '<div class="col-xs-3">' + '<div class="table__th">Category / Node / Property</div>' + '</div>' + '<div class="col-xs-9">' + '<div class="table__thead row">' + '<div class="table__th col-xs-6">Matched GDC Values <a class="table__tooltip tooltip-target" data-toggle="tooltip" data-placement="bottom" title="Values that are found in the GDC dictionary and may be successfully submitted for the corresponding property."><i class="fa fa-info-circle"></i></a></div>' + '<div class="table__th col-xs-6">CDE Permissible Values <a class="table__tooltip tooltip-target" data-toggle="tooltip" data-placement="bottom" title="For GDC dictionary properties that have a corresponding caDSR clinical data element (CDE), these values are part of that CDE\'s value domain in the caDSR. They may not currently be available in the GDC dictionary."><i class="fa fa-info-circle"></i></a></div>' + '</div>' + '</div>' + '</div>' + '<div id="table-body" class="table__body table__body--overflow row" style="max-height: {{:mh}}px;"><div class="col-xs-12">{{for values}}' + '<div class="table__row row table__row--striped table__row--flex">' + '<div class="table__td col-xs-3">' + '{{:category}}<ul class="table__ul">' + '<li class="table__li table__td--word-break">{{:node}}' + '<ul class="table__ul"><li class="table__li table__td--word-break">{{:name}}</li></ul>' + '</li></ul>' + '<a href="javascript:void(0)" class="gdc-details"><i class="fa fa-angle-down"></i> detail</a>' + '<div id="gdc-links" style="display: none;">' + '{{if local}}' + '<a href="javascript:getGDCData(\'{{:ref}}\',null);">See All Values</a></br>' + '<a href="javascript:toCompare(\'{{:ref}}\');"> Compare with User List</a></br>' + '{{/if}}' + '{{if syn}}' + '<a href="javascript:getGDCTerms(\'{{:ref}}\', \'{{:tgts_enum_n}}\');">See All Terms</a></br>' + '{{/if}}' + '{{if cdeId == ""}}' + '' + '{{else}}' + 'caDSR: <a class="table-td-link" href="{{:cdeUrl}}" target="_blank">CDE</a>' + '{{if local && cdeLen}}' + ' , <a class="table-td-link" href="javascript:getCDEData(\'{{:cdeId}}\', \'{{:tgts_cde_n}}\');">Values</a> , <a class="table-td-link" href="javascript:compareGDC(\'{{:ref}}\',\'{{:cdeId}}\');"> Compare with GDC</a>' + '{{else cdeLen}}' + ' , <a class="table-td-link" href="javascript:getCDEData(\'{{:cdeId}}\', \'{{:tgts_cde_n}}\');">Values</a>' + '{{else}}' + '' + '{{/if}}' + '{{/if}}' + '</div>' + '</div>' + '<div class="table__values col-xs-9"> {{for vs}}' + '<div class="row table__row--flex{{if #getIndex() >= 5}} table__row--toggle{{/if}}">' + '<div class="table__td table__gdc-values col-xs-6">' + '{{if n == "no match"}}' + '<div class="row">' + '<div class="col-xs-9">no match</div>' + '<div class="col-xs-3"><a href="javascript:void(0);" class="cde-suggest" style="float: right;">Suggest Item</a></div>' + '</div>' + '{{else}}' + '{{if i_c !== undefined }}' + '<div class="row">' + '<div class="col-xs-10">' + '<div class="row">' + '<div class="col-xs-2 table__ico3-code">' + '{{:i_c}}' + '</div>' + '<div class="col-xs-10">' + '{{for n ~ref=ref}}' + '<a href="javascript:getGDCData(\'{{:~ref}}\',\'{{:}}\');">{{:}} (ICD-O-3)</a><br>' + '{{/for}}' + '</div>' + '</div>' + '</div>' + '<div class="col-xs-2 table__collapser">{{if n_t.length }}<a href="javascript:void(0);" class="collapser" aria-label="collapser"><i class="fa fa-plus"></i></a>{{/if}}</div>' + '</div>' + '<div id="data-content" class="table__td" style="display: none;">' + '{{for n_t}}' + '<div class="row table__row-syn">' + '<div class="col-xs-4">' + '<a class="table-td-link" href="javascript:getNCITDetails(\'{{:n_c}}\');">{{:n_c}}</a> (NCIt)' + '</div>' + '<div class="col-xs-8">{{for s}}{{:}}</br>{{/for}}</div>' + '</div>' + '{{/for}}' + '</div>' + '{{else}}' + '<div class="row">' + '<div class="col-xs-10"><a href="javascript:getGDCData(\'{{:ref}}\',\'{{:n}}\');">' + '{{if i_c !== undefined }}{{:i_c}} {{:n}} (ICD-O-3){{else}}{{:n}}{{/if}}' + '</a></div>' + '<div class="col-xs-2 table__collapser">{{if s.length }}<a href="javascript:void(0);" class="collapser" aria-label="collapser"><i class="fa fa-plus"></i></a>{{/if}}</div>' + '</div>' + '<div id="data-content" class="table__td" style="display: none;">' + '<div class="row">' + '<div class="col-xs-4">' + '{{* if((/^C[1-9]/g).test(data.n_c)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:n_c}}\');">{{:n_c}}</a> (NCIt)' + '{{* } else { }} {{:c}} {{*: (/^[E]/g).test(data.n_c) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '</div>' + '<div class="col-xs-8">{{for s}}{{:}}</br>{{/for}}</div>' + '</div>' + '</div>' + '{{/if}}' + '{{/if}}' + '</div>' + '<div class="table__td table__cde-values col-xs-6">' + '{{if cde_s.length }}' + '<div class="row">' + '<div class="col-xs-10">{{:cde_pv}}</div>' + '<div class="col-xs-2 table__collapser">' + '<a href="javascript:void(0);" class="collapser" aria-label="collapser"><i class="fa fa-plus"></i></a>' + '</div>' + '</div>' + '<div id="data-content" class="table__td" style="display: none;">' + '<div class="row">' + '<div class="table__td col-xs-12">PV Meaning (caDSR): {{:cde_pvm}}</div>' + '</div>' + '{{for cde_s}}' + '<div class="row">' + '<div class="col-xs-4">' + '{{* if((/^C[1-9]/g).test(data.c)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:c}}\');">{{:c}}</a> (NCIt)' + '{{* } else { }} {{:c}} {{*: (/^[E]/g).test(data.c) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '</div>' + '<div class="col-xs-8">{{for s}}{{:}}</br>{{/for}}</div>' + '</div>' + '{{/for}}' + '</div>' + '{{/if}}' + '</div>' + '</div> {{/for}}' + '{{if vs.length > 5}}' + '<div class="row"><div class="table__td col-xs-12">' + '<a class="table-td-link show-more-less" href="javascript:void(0);"><i class="fa fa-angle-down"></i> Show More ({{:vs.length - 5}})</a>' + '</div></div>' + '{{/if}}' + '</div>' + '</div> {{/for}} </div></div></div>' + '<div id="alert-suggest" class="alert alert__suggest alert-info alert-dismissible" role="alert" style="display: none;">' + 'An email will be sucessfully sent to <strong>GDC</strong> and <strong>EVS</strong> team.' + '</div>';

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 12 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(13);


/* harmony default export */ __webpack_exports__["a"] = (function (trsHtml, psHtml, vsHtml, keyword) {

  var html = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */]).render({
    trs_active: trsHtml.active,
    trs_len: trsHtml.len,
    trsHtml: trsHtml.html,
    ps_active: psHtml.active,
    ps_len: psHtml.len,
    psHtml: psHtml.html,
    vs_active: vsHtml.active,
    vs_len: vsHtml.len,
    vsHtml: vsHtml.html,
    keyword: keyword
  });

  return html;
});

/***/ }),
/* 13 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

var tmpl = '<div class="results">' + '<div class="tab-nav">' + '<div class="tab-nav__text">Results for <span class="tab-nav__term">\'{{:keyword}}\'</span> in:</div>' + '<div id="tab" class="btn-group tab-nav__group" data-toggle="buttons">' + '<ul class="tab-nav__ul"role="tablist">' + '<li role="presentation" class="tab-nav__li{{if vs_active}} active{{else}}{{/if}}">' + '<div class="tab-nav__tooltip tooltip-target" data-toggle="tooltip" data-placement="bottom" data-delay=\'{"show":"1000"}\' data-trigger="hover" title="Text will be provided to inform users on how to interpret content of tabs.">' + '<a id="tab-values" href="#values" class="tab-nav__btn" aria-controls="values" role="tab" data-toggle="tab" aria-expanded="true">Values</a>' + '<span class="tab-nav__notification">{{:vs_len}}</span>' + '</div>' + '</li>' + '<li role="presentation" class="tab-nav__li{{if ps_active}} active{{else}}{{/if}}">' + '<div class="tab-nav__tooltip tooltip-target" data-toggle="tooltip" data-placement="bottom" data-delay=\'{"show":"1000"}\' data-trigger="hover" title="Text will be provided to inform users on how to interpret content of tabs.">' + '<a id="tab-properties" href="#properties" class="tab-nav__btn" aria-controls="properties" role="tab" data-toggle="tab" aria-expanded="true">Properties</a>' + '<span class="tab-nav__notification">{{:ps_len}}</span>' + '</div>' + '</li>' + '<li role="presentation" class="tab-nav__li{{if trs_active}} active{{else}}{{/if}}">' + '<div class="tab-nav__tooltip tooltip-target" data-toggle="tooltip" data-placement="bottom" data-delay=\'{"show":"1000"}\' data-trigger="hover" title="Text will be provided to inform users on how to interpret content of tabs.">' + '<a id="tab-dictionary" href="#dictionary" class="tab-nav__btn" aria-controls="dictionary" role="tab" data-toggle="tab" aria-expanded="true">Dictionary</a>' + '<span class="tab-nav__notification">{{:trs_len}}</span>' + '</div>' + '</li>' + '</ul>' + '</div></div>' + '<div class="tab-content"><div role="tabpanel" class="tab-pane {{if vs_active}}active{{else}}{{/if}}" id="values">{{:vsHtml}}</div>' + '<div role="tabpanel" class="tab-pane {{if ps_active}}active{{else}}{{/if}}" id="properties">{{:psHtml}}</div>' + '<div role="tabpanel" class="tab-pane {{if trs_active}}active{{else}}{{/if}}" id="dictionary">{{:trsHtml}}</div></div></div>';

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 14 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

var tmpl = '{{for options}}<div class="suggest__object">' + '<span class="suggest__name">{{:id}}</span>' + '<label class="suggest__type">{{:type}}</label>' + '</div>{{/for}}';

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 15 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(16);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__api__ = __webpack_require__(0);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__shared__ = __webpack_require__(1);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_3__gdc_data__ = __webpack_require__(17);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_4__gdc_terms__ = __webpack_require__(19);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_5__cde_data__ = __webpack_require__(21);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_6__ncit_details__ = __webpack_require__(23);








var func = {
    getGDCData: function getGDCData(prop, item) {
        Object(__WEBPACK_IMPORTED_MODULE_3__gdc_data__["a" /* default */])(prop, item);
    },
    getGDCTerms: function getGDCTerms(uid, tgts) {
        Object(__WEBPACK_IMPORTED_MODULE_4__gdc_terms__["a" /* default */])(uid, tgts);
    },
    toCompare: function toCompare(uid) {
        __WEBPACK_IMPORTED_MODULE_1__api__["a" /* default */].getGDCDataById(uid, function (id, items) {
            if ($('#compare_dialog').length) {
                $('#compare_dialog').remove();
            }
            var windowEl = $(window);
            var html = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].toCompare).render({ items: items });

            var tp = window.innerHeight * 0.2 < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() ? 20 : window.innerHeight * 0.2;
            //display result in a table
            $(document.body).append(html);
            $("#compare_dialog").dialog({
                modal: false,
                position: { my: "center top+" + tp, at: "center top", of: $('#docs-container') },
                //width:"60%",
                width: 750,
                height: 630,
                minWidth: 750,
                maxWidth: 900,
                minHeight: 542,
                maxHeight: 800,
                title: "Compare Your Values with GDC Values ",
                open: function open() {

                    var target = $(this).parent();
                    target.find('.ui-dialog-titlebar').css('padding', '15px');
                    target.find('.ui-dialog-titlebar-close').html('');
                    if (target.offset().top - windowEl.scrollTop() < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset()) {
                        target.css('top', windowEl.scrollTop() + __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 + 'px');
                    }

                    $('#cp_result').css("display", "none");
                    $('#compare').bind('click', function () {
                        var gv = [];
                        items.forEach(function (item) {
                            gv.push(item.n);
                        });
                        compare(gv);
                    });
                    $('#cancelCompare').bind('click', function () {
                        $("#compare_dialog").dialog('close');
                    });
                },
                close: function close() {
                    $(this).remove();
                }
            }).parent().draggable({
                containment: '#docs-container'
            });
        });
    },
    getCDEData: function getCDEData(uid, tgts) {
        Object(__WEBPACK_IMPORTED_MODULE_5__cde_data__["a" /* default */])(uid, tgts);
    },
    compareGDC: function compareGDC(prop, uid) {
        var ids = {};
        ids.local = prop;
        ids.cde = uid;
        __WEBPACK_IMPORTED_MODULE_1__api__["a" /* default */].getGDCandCDEDataById(ids, function (ids, items) {
            if ($('#compareGDC_dialog').length) {
                $('#compareGDC_dialog').remove();
            }
            var windowEl = $(window);
            var popup = '<div id="compareGDC_dialog">' + '<div id="compareGDC_result"></div>' + '</div>';
            $(document.body).append(popup);
            var tp = window.innerHeight * 0.2 < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() ? __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 : window.innerHeight * 0.2;
            var toV = [];
            var fromV = [];
            var opt = {};
            opt.sensitive = false;
            opt.unmatched = false;
            items.to.forEach(function (t) {
                toV.push(t.n);
            });
            items.from.forEach(function (f) {
                fromV.push(f.n);
            });
            var table = generateCompareGDCResult(fromV, toV, opt);
            var html = '<div id="cpGDC_result_option">'
            //+'<div class="option-left"><input type="checkbox" id="compareGDC_filter"> Case Sensitive</div><div class="option-right"><input type="checkbox" id="compareGDC_unmatched"> Hide Unmatched Values</div></div><div class="clearfix"></div>'
            + '<div id="cpGDC_result_table" class="table__container">' + table + '</div>' + '</div>';

            var titleComponent = '<div class="checkbox ui-checkbox"><label class="checkbox__label checkbox__label--height"><input id="compareGDC_filter" class="checkbox__input" type="checkbox" value=""><span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Case Sensitive</label>' + '<label class="checkbox__label checkbox__label--height"><input id="compareGDC_unmatched" class="checkbox__input" type="checkbox" value=""><span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Hide Unmatched Values</label></div>';

            $('#compareGDC_result').html(html);

            $("#compareGDC_dialog").dialog({
                modal: false,
                position: { my: "center top+" + tp, at: "center top", of: $('#docs-container') },
                //width:"50%",
                width: 750,
                height: 550,
                minWidth: 715,
                maxWidth: 900,
                minHeight: 300,
                maxHeight: 800,
                title: "Compare GDC Values with caDSR Values ",
                open: function open() {

                    var target = $(this).parent();
                    target.find('.ui-dialog-titlebar').css('padding', '15px').append(titleComponent);
                    target.find('.ui-dialog-titlebar-close').html('');
                    if (target.offset().top - windowEl.scrollTop() < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset()) {
                        target.css('top', windowEl.scrollTop() + __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 + 'px');
                    }

                    $('#compareGDC_filter').bind('click', function () {
                        var options = {};
                        options.sensitive = $("#compareGDC_filter").prop('checked');
                        options.unmatched = $("#compareGDC_unmatched").prop('checked');
                        var table_new = generateCompareGDCResult(fromV, toV, options);
                        $('#cpGDC_result_table').html(table_new);
                    });
                    $('#compareGDC_unmatched').bind('click', function () {
                        var options = {};
                        options.sensitive = $("#compareGDC_filter").prop('checked');
                        options.unmatched = $("#compareGDC_unmatched").prop('checked');
                        var table_new = generateCompareGDCResult(fromV, toV, options);
                        $('#cpGDC_result_table').html(table_new);
                    });
                },
                close: function close() {
                    $(this).remove();
                }
            }).parent().draggable({
                containment: '#docs-container'
            });
        });
    },
    getNCITDetails: function getNCITDetails(uid) {
        Object(__WEBPACK_IMPORTED_MODULE_6__ncit_details__["a" /* default */])(uid);
    }
};

/* harmony default export */ __webpack_exports__["a"] = (func);

/***/ }),
/* 16 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

var tmpl = {
    toCompare: '<div id="compare_dialog">' + '<div id="compare_form" class="compare-form">' + '<div id="cp_top" class="compare-form__top">' + '<label class="compare-form__label--left">User Defined Values:</label>' + '<label class="compare-form__label--right">GDC Values:</label>' + '<div id="cp_left" class="compare-form__left">' + '<textarea id="cp_input" class="compare-form__textarea" rows="10" cols="20" placeholder="Input values line by line" autocomplete="off"></textarea></div>' + '<div id="cp_middle" class="compare-form__middle"></div>' + '<div id="cp_right" class="compare-form__right">' + '{{for items}}' + '<div>{{:n}}</div>' + '{{/for}}' + '</div>' + '</div>' + '<div id="cp_massage" class="compare-form__message"></div>' + '<div id="cp_bottom" class="compare-form__bottom">' + '<button id="compare" class="btn btn-default compare-form__button">Compare</button>' + '<button id="cancelCompare" class="btn btn-default compare-form__button">Cancel</button>' + '</div>' + '</div>' + '<div id="compare_result" class="compare_result"></div>' + '</div>'
};

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 17 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony export (immutable) */ __webpack_exports__["a"] = gdcData;
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(18);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__api__ = __webpack_require__(0);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__shared__ = __webpack_require__(1);




function gdcData(prop, item) {
  __WEBPACK_IMPORTED_MODULE_1__api__["a" /* default */].getGDCDataById(prop, function (id, items) {

    if ($('#gdc_data').length) {
      $('#gdc_data').remove();
    }

    var windowEl = $(window);
    var icdo = false;
    items.forEach(function (item) {
      if (item.i_c !== undefined) {
        icdo = true;
      }
    });

    var target = item == undefined ? item : item.replace(/<b>/g, "").replace(/<\/b>/g, "");
    var header = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].header).render({ target: target, icdo: icdo, items_length: items.length });
    var html = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].body).render({ target: target, icdo: icdo, items: items });
    var tp = window.innerHeight * 0.2 < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() ? 20 : window.innerHeight * 0.2;

    //display result in a table
    $(document.body).append(html);

    $('#gdc_data').dialog({
      modal: false,
      position: { my: 'center top+' + tp, at: 'center top', of: $('#docs-container') },
      width: 600,
      height: 450,
      minWidth: 420,
      maxWidth: 800,
      minHeight: 350,
      maxHeight: 650,
      open: function open() {
        //add new custom header
        if (icdo) {
          $(this).prev('.ui-dialog-titlebar').css('padding-top', '7.5em').html(header);
        } else {
          $(this).prev('.ui-dialog-titlebar').css('padding-top', '3.5em').html(header);
        }

        var target = $(this).parent();
        if (target.offset().top - windowEl.scrollTop() < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset()) {
          target.css('top', windowEl.scrollTop() + __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 + 'px');
        }

        $('#close_gdc_data').bind('click', function () {
          $("#gdc_data").dialog('close');
        });
      },
      close: function close() {
        $(this).remove();
      }
    }).parent().draggable({
      containment: '#docs-container'
    });

    if ($('#show_all_gdc_data') !== undefined) {
      $('#show_all_gdc_data').bind('click', function () {
        var v = $(this).prop("checked");
        if (v) {
          $('#gdc-data-list div[style="display: none;"]').each(function () {
            $(this).css("display", "block");
          });
          var setScroll = $('#gdc_data_match').offset().top - $('#gdc_data').offset().top;
          $('#gdc_data').scrollTop(setScroll - 120);
        } else {
          $('#gdc-data-list div[style="display: block;"]').each(function () {
            $(this).css("display", "none");
          });
        }
      });
    }
  });
}

/***/ }),
/* 18 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
var tmpl = {
  header: '<div class="dialog__header">' + '<div class="dialog__titlebar">' + '<span id="ui-id-4" class="ui-dialog-title">GDC Values</span>' + '<button type="button" id="close_gdc_data" class="ui-button ui-corner-all ui-widget ui-button-icon-only ui-dialog-titlebar-close" title="Close"></button>' + '<span class="ui-label">{{:items_length}}</span>' + '{{if target !== null}}' + '<div class="checkbox ui-checkbox">' + '<label class="checkbox__label checkbox__label--height">' + '<input id="show_all_gdc_data" class="checkbox__input" type="checkbox" value="">' + '<span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Show all GDC values' + '</label>' + '</div>' + '{{/if}}' + '</div>' + '{{if icdo !== false}}' + '<div class="table__container">' + '<div class="table__thead row">' + '<div class="table__th col-xs-4">ICD-O-3 Code</div>' + '<div class="table__th col-xs-8">ICD-O-3 Term</div>' + '</div>' + '</div>' + '{{/if}}' + '</div>',
  body: '<div id="gdc_data">' + '{{if icdo}}' + '<div class="table__container">' + '<div id="gdc-data-list" class="table__body row">' + '<div class="col-xs-12">' + '{{if target !== null }}' + '{{for items}}' + '{{if n == ~root.target }}' + '<div class="table__row row" id="gdc_data_match">' + '<div class="table__td table__td--slim col-xs-4">{{:i_c.c}}</div><div class="table__td table__td--slim col-xs-8"><b>{{:n}}</b></div>' + '</div>' + '{{else}}' + '<div class="table__row row" style="display: none;">' + '<div class="table__td table__td--slim col-xs-4">{{:i_c.c}}</div><div class="table__td table__td--slim col-xs-8">{{:n}}</div>' + '</div>' + '{{/if}}' + '{{/for}}' + '{{else}}' + '{{for items}}' + '<div class="table__row row">' + '<div class="table__td table__td--slim col-xs-4">{{:i_c.c}}</div><div class="table__td table__td--slim col-xs-8">{{:n}}</div>' + '</div>' + '{{/for}}' + '{{/if}}' + '</div>' + '</div>' + '</div>' + '{{else}}' + '<div class=" table__container table__container--blank">' + '<div id="gdc-data-list" class="table__body row">' + '<div class="col-xs-12">' + '{{if target !== null }}' + '{{for items}}' + '{{if n == ~root.target }}' + '<div class="row" id="gdc_data_match">' + '<div class="table__td table__td--xslim col-xs-8"><b>{{:n}}</b></div>' + '</div>' + '{{else}}' + '<div class="row" style="display: none;">' + '<div class="table__td table__td--xslim col-xs-8">{{:n}}</div>' + '</div>' + '{{/if}}' + '{{/for}}' + '{{else}}' + '{{for items}}' + '<div class="row">' + '<div class="table__td table__td--xslim col-xs-8">{{:n}}</div>' + '</div>' + '{{/for}}' + '{{/if}}' + '</div>' + '</div>' + '</div>' + '{{/if}}' + '</div>'
};

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 19 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony export (immutable) */ __webpack_exports__["a"] = getGDCTerms;
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(20);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__api__ = __webpack_require__(0);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__shared__ = __webpack_require__(1);




function getGDCTerms(uid, tgts) {
  __WEBPACK_IMPORTED_MODULE_1__api__["a" /* default */].getGDCDataById(uid, function (id, items) {
    if ($('#gdc_terms_data').length) {
      $('#gdc_terms_data').remove();
    }

    var targets = null;
    var icdo = false;
    var windowEl = $(window);

    if (tgts !== null && tgts !== undefined) {
      targets = tgts.split("#");

      items.forEach(function (item) {
        if (item.i_c !== undefined) {
          icdo = true;
        }
        if (targets.indexOf(item.n) > -1) {
          item.e = true;
        }
      });
    } else {
      items.forEach(function (item) {
        if (item.i_c !== undefined) {
          icdo = true;
        }
      });
    }

    items.forEach(function (it) {
      if (it.s == undefined) return;
      var cache = {};
      var tmp_s = [];
      it.s.forEach(function (s) {
        var lc = s.trim().toLowerCase();
        if (!(lc in cache)) {
          cache[lc] = [];
        }
        cache[lc].push(s);
      });
      for (var idx in cache) {
        //find the term with the first character capitalized
        var word = findWord(cache[idx]);
        tmp_s.push(word);
      }
      it.s_r = tmp_s;
    });

    var header = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].header).render({ targets: targets, icdo: icdo, items_length: items.length });
    var html = $.templates({ markup: __WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].body, allowCode: true }).render({ targets: targets, icdo: icdo, items: items });
    var tp = window.innerHeight * 0.2 < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() ? __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 : window.innerHeight * 0.2;

    //display result in a table
    $(document.body).append(html);

    $("#gdc_terms_data").dialog({
      modal: false,
      position: { my: "center top+" + tp, at: "center top", of: $('#docs-container') },
      width: 900,
      height: 'auto',
      minWidth: 700,
      maxWidth: 1000,
      minHeight: 300,
      maxHeight: 600,
      open: function open() {
        //add new custom header
        $(this).prev('.ui-dialog-titlebar').css('padding-top', '7.5em').html(header);

        // $(this).prev('.ui-dialog-titlebar').remove();
        // $(this).before(header);
        //$(this).before(header);
        var target = $(this).parent();

        // target.find('.ui-dialog-titlebar').append(titleComponent);
        // target.find('.ui-dialog-titlebar-close').html('');

        if (target.offset().top - windowEl.scrollTop() < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset()) {
          target.css('top', windowEl.scrollTop() + __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 + 'px');
        }

        $('#close_gdc_terms_data').bind('click', function () {
          $("#gdc_terms_data").dialog('close');
        });

        $('#gdc-data-invariant').bind('click', function () {
          $("#gdc-syn-data-list").find('div[name="syn_area"]').each(function () {
            var rp = $(this).html();
            var invariant = $(this).parent().children('div[name="syn_invariant"]');
            $(this).html(invariant[0].innerHTML);
            invariant[0].innerHTML = rp;
          });
        });
      },
      close: function close() {
        $(this).remove();
      }
    }).parent().draggable({
      containment: '#docs-container'
    });

    if ($('#show_all_gdc_syn') !== undefined) {
      $('#show_all_gdc_syn').bind('click', function () {
        var v = $(this).prop("checked");
        if (v) {
          $('#gdc-syn-data-list div.table__row[style="display: none;"]').each(function () {
            $(this).css("display", "block");
          });
        } else {
          $('#gdc-syn-data-list div.table__row[style="display: block;"]').each(function () {
            $(this).css("display", "none");
          });
        }
      });
    }
  });
}

/***/ }),
/* 20 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
var tmpl = {
  header: '<div class="dialog__header">' + '<div class="dialog__titlebar">' + '<span id="ui-id-4" class="ui-dialog-title">GDC Terms</span>' + '<button type="button" id="close_gdc_terms_data" class="ui-button ui-corner-all ui-widget ui-button-icon-only ui-dialog-titlebar-close" title="Close"></button>' + '<span class="ui-label">{{:items_length}}</span>' + '<div class="checkbox ui-checkbox">' + '<label class="checkbox__label checkbox__label--height">' + '<input id="gdc-data-invariant" class="checkbox__input" type="checkbox" value="">' + '<span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Show Duplicates' + '</label>' + '{{if targets !== null}}' + '<label class="checkbox__label checkbox__label--height">' + '<input id="show_all_gdc_syn" class="checkbox__input" type="checkbox" value="">' + '<span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Show all GDC values' + '</label>' + '{{/if}}' + '</div>' + '</div>' + '<div class="table__container">' + '<div class="table__thead row">' + '{{if icdo}}' + '<div class="table__th col-xs-2">ICD-O-3 Code</div>' + '<div class="table__th col-xs-3">ICD-O-3 Term</div>' + '{{else}}' + '<div class="table__th col-xs-5">GDC Term</div>' + '{{/if}}' + '<div class="table__th col-xs-2">NCIt Code</div>' + '<div class="table__th col-xs-5">NCIt Terms</div>' + '</div>' + '</div>' + '</div>',
  body: '<div id="gdc_terms_data">' + '<div id="gdc-syn-data-list" class="table__container">' + '<div class="table__body row">' + '<div class="col-xs-12">' + '{{for items}}' + '{{if e == true || ~root.targets == null}}' + '<div class="table__row row">' + '{{if ~root.icdo}}' + '<div class="table__td col-xs-2">{{:i_c.c}}</div>' + '<div class="table__td col-xs-3">{{:n}}</div>' + '{{else}}' + '<div class="table__td col-xs-5">{{:n}}</div>' + '{{/if}}' + '<div class="table__td col-xs-2">' + '{{if n_c && n_c !== ""}}' + '{{* if((/^C[1-9]/g).test(data.n_c)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:n_c}}\');">{{:n_c}}</a> (NCIt)' + '{{* } else { }} {{:n_c}} {{*: (/^[E]/g).test(data.n_c) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '{{/if}}' + '</div>' + '<div name="syn_area" class="table__td col-xs-5">{{for s_r}}{{>#data}}<br>{{/for}}</div>' + '<div name="syn_invariant" class="table__td col-xs-5" style="display: none;">' + '{{for s}}{{>#data}}<br>{{/for}}' + '</div>' + '</div>' + '{{else}}' + '<div class="table__row row" style="display: none;">' + '{{if ~root.icdo}}' + '<div class="table__td col-xs-2">{{:i_c.c}}</div>' + '<div class="table__td col-xs-3">{{:n}}</div>' + '{{else}}' + '<div class="table__td col-xs-5">{{:n}}</div>' + '{{/if}}' + '<div class="table__td col-xs-2">' + '{{if n_c && n_c !== ""}}' + '{{* if((/^C[1-9]/g).test(data.n_c)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:n_c}}\');">{{:n_c}}</a> (NCIt)' + '{{* } else { }} {{:n_c}} {{*: (/^[E]/g).test(data.n_c) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '{{/if}}' + '</div>' + '<div name="syn_area" class="table__td col-xs-5">{{for s_r}}{{>#data}}<br>{{/for}}</div>' + '<div name="syn_invariant" class="table__td col-xs-5" style="display: none;">' + '{{for s}}{{>#data}}<br>{{/for}}' + '</div>' + '</div>' + '{{/if}}' + '{{/for}}' + '</div>' + '</div>' + '</div>' + '</div>'
};

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 21 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony export (immutable) */ __webpack_exports__["a"] = cdeData;
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(22);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__api__ = __webpack_require__(0);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__shared__ = __webpack_require__(1);




function cdeData(uid, tgts) {
  __WEBPACK_IMPORTED_MODULE_1__api__["a" /* default */].getCDEDataById(uid, function (id, items) {
    //data precessing
    var tmp = [];
    items.forEach(function (item) {
      var t = {};
      t.pv = item.n;
      t.pvm = item.m;
      t.pvd = item.d;
      t.i_rows = [];
      t.rows = [];
      item.ss.forEach(function (s) {
        var i_r = {};
        var r = {};
        i_r.pvc = s.c;
        r.pvc = s.c;
        r.s = s.s;
        i_r.s = [];
        //remove duplicate
        var cache = {};
        s.s.forEach(function (w) {
          var lc = w.trim().toLowerCase();
          if (!(lc in cache)) {
            cache[lc] = [];
          }
          cache[lc].push(w);
        });
        for (var idx in cache) {
          //find the term with the first character capitalized
          var word = findWord(cache[idx]);
          i_r.s.push(word);
        }
        t.i_rows.push(i_r);
        t.rows.push(r);
      });
      tmp.push(t);
    });

    var targets = null;
    var windowEl = $(window);

    if (tgts !== null && tgts !== undefined && tgts !== "") {
      tgts = tgts.replace(/\^/g, '\'');
      targets = tgts.split("#");

      tmp.forEach(function (item) {
        if (targets.indexOf(item.pv) > -1) {
          item.e = true;
        }
      });
    }

    if ($('#caDSR_data').length) {
      $('#caDSR_data').remove();
    }

    var header = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].header).render({ targets: targets, items_length: items.length });

    var html = $.templates({ markup: __WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].body, allowCode: true }).render({ targets: targets, items: tmp });
    var tp = window.innerHeight * 0.2 < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() ? __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 : window.innerHeight * 0.2;

    //display result in a table
    $(document.body).append(html);

    $("#caDSR_data").dialog({
      modal: false,
      position: { my: "center top+" + tp, at: "center top", of: $('#docs-container') },
      width: 900,
      height: 'auto',
      minWidth: 700,
      maxWidth: 1000,
      minHeight: 300,
      maxHeight: 600,
      title: "caDSR Values",
      open: function open() {
        //add new custom header
        $(this).prev('.ui-dialog-titlebar').css('padding-top', '7.5em').html(header);

        var target = $(this).parent();
        if (target.offset().top - windowEl.scrollTop() < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset()) {
          target.css('top', windowEl.scrollTop() + __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 + 'px');
        }

        $('#close_caDSR_data').bind('click', function () {
          $("#caDSR_data").dialog('close');
        });

        $('#cde-data-invariant').bind('click', function () {
          $("#cde-syn-data-list").find('div[name="syn_area"]').each(function () {
            var rp = $(this).html();
            var invariant = $(this).parent().children('div[name="syn_invariant"]');
            $(this).html(invariant[0].innerHTML);
            invariant[0].innerHTML = rp;
          });
        });
      },
      close: function close() {
        $(this).remove();
      }
    }).parent().draggable({
      containment: '#docs-container'
    });

    if ($('#show_all_cde_syn') !== undefined) {
      $('#show_all_cde_syn').bind('click', function () {
        var v = $(this).prop("checked");
        if (v) {
          $('#cde-syn-data-list div.table__row[style="display: none;"]').each(function () {
            $(this).css("display", "block");
          });
        } else {
          $('#cde-syn-data-list div.table__row[style="display: block;"]').each(function () {
            $(this).css("display", "none");
          });
        }
      });
    }
  });
}

/***/ }),
/* 22 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
var tmpl = {
  header: '<div class="dialog__header">' + '<div class="dialog__titlebar">' + '<span id="ui-id-4" class="ui-dialog-title">caDSR Values</span>' + '<button type="button" id="close_caDSR_data" class="ui-button ui-corner-all ui-widget ui-button-icon-only ui-dialog-titlebar-close" title="Close"></button>' + '<span class="ui-label">{{:items_length}}</span>' + '<div class="checkbox ui-checkbox">' + '<label class="checkbox__label checkbox__label--height">' + '<input id="cde-data-invariant" class="checkbox__input" type="checkbox" value="">' + '<span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Show Duplicates' + '</label>' + '{{if targets !== null}}' + '<label class="checkbox__label checkbox__label--height">' + '<input id="show_all_cde_syn" class="checkbox__input" type="checkbox" value="">' + '<span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Show all CDE values' + '</label>' + '{{/if}}' + '</div>' + '</div>' + '<div class="table__container">' + '<div class="table__thead row">' + '<div class="table__th col-xs-2">PV</div>' + '<div class="table__th col-xs-2">PV Meaning</div>' + '<div class="table__th col-xs-4">Description</div>' + '<div class="table__th col-xs-4">NCIt Code and Synonyms</div>' + '</div>' + '</div>' + '</div>',
  body: '<div id="caDSR_data">' + '<div id="cde-syn-data-list" class="table__container">' + '<div class="table__body row">' + '<div class="col-xs-12">' + '{{for items}}' + '{{if e == true || ~root.targets == null}}' + '<div class="table__row row">' + '<div class="table__td col-xs-2">{{:pv}}</div>' + '<div class="table__td col-xs-2">{{:pvm}}</div>' + '<div class="table__td col-xs-4">{{:pvd}}</div>' + '<div name="syn_area" class="table__td col-xs-4">' + '{{for i_rows}}' + '<div class="row">' + '<div class="col-lg-3">' + '{{* if((/^C[1-9]/g).test(data.pvc)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:pvc}}\');">{{:pvc}}</a> (NCIt)' + '{{* } else { }} {{:pvc}} {{*: (/^[E]/g).test(data.pvc) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '</div>' + '<div class="col-lg-9 col-xs-12">{{for s}}{{>#data}}<br>{{/for}}</div>' + '</div>' + '{{/for}}' + '</div>' + '<div name="syn_invariant" class="table__td col-xs-4" style="display: none;">' + '{{for rows}}' + '<div class="row">' + '<div class="col-lg-3">' + '{{* if((/^C[1-9]/g).test(data.pvc)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:pvc}}\');">{{:pvc}}</a> (NCIt)' + '{{* } else { }} {{:pvc}} {{*: (/^[E]/g).test(data.pvc) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '</div>' + '<div class="col-lg-9 col-xs-12">{{for s}}{{>#data}}<br>{{/for}}</div>' + '</div>' + '{{/for}}' + '</div>' + '</div>' + '{{else}}' + '<div class="table__row row" style="display: none;">' + '<div class="table__td col-xs-2">{{:pv}}</div>' + '<div class="table__td col-xs-2">{{:pvm}}</div>' + '<div class="table__td col-xs-4">{{:pvd}}</div>' + '<div name="syn_area" class="table__td col-xs-4">' + '{{for i_rows}}' + '<div class="row">' + '<div class="col-lg-3 col-xs-12">' + '{{* if((/^C[1-9]/g).test(data.pvc)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:pvc}}\');">{{:pvc}}</a> (NCIt)' + '{{* } else { }} {{:pvc}} {{*: (/^[E]/g).test(data.pvc) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '</div>' + '<div class="col-lg-9 col-xs-12">{{for s}}{{>#data}}<br>{{/for}}</div>' + '</div>' + '{{/for}}' + '</div>' + '<div name="syn_invariant" class="table__td col-xs-4" style="display: none;">' + '{{for rows}}' + '<div class="row">' + '<div class="col-lg-3 col-xs-12">' + '{{* if((/^C[1-9]/g).test(data.pvc)) { }}<a class="table-td-link" href="javascript:getNCITDetails(\'{{:pvc}}\');">{{:pvc}}</a> (NCIt)' + '{{* } else { }} {{:pvc}} {{*: (/^[E]/g).test(data.pvc) ? "(CTCAE)" : "(NCIt)" }} {{* } }}' + '</div>' + '<div class="col-lg-9 col-xs-12">{{for s}}{{>#data}}<br>{{/for}}</div>' + '</div>' + '{{/for}}' + '</div>' + '</div>' + '{{/if}}' + '{{/for}}' + '</div>' + '</div>' + '</div>' + '</div>',
  footer: ''
};

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 23 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony export (immutable) */ __webpack_exports__["a"] = ncitDetails;
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__view__ = __webpack_require__(24);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__api__ = __webpack_require__(0);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__shared__ = __webpack_require__(1);




function ncitDetails(uid) {
  __WEBPACK_IMPORTED_MODULE_1__api__["a" /* default */].evsRestApi(uid, function (id, item) {

    if ($('#ncit_details').length) {
      $('#ncit_details').remove();
    }

    var tmp = {};
    tmp.code = item.code;
    tmp.name = item.preferredName;
    tmp.definition = item.definitions.length ? item.definitions.find(function (defs) {
      return defs.defSource === 'NCI';
    }).description : undefined;
    var tmp_s = item.synonyms.map(function (syns) {
      return syns.termName;
    });
    tmp.synonyms = [];
    //remove the duplicate

    var cache = {};
    if (tmp_s.length > 0) {
      tmp_s.forEach(function (s) {
        var lc = s.trim().toLowerCase();
        if (!(lc in cache)) {
          cache[lc] = [];
        }
        cache[lc].push(s);
      });
      for (var idx in cache) {
        //find the term with the first character capitalized
        var word = findWord(cache[idx]);
        tmp.synonyms.push(word);
      }
    }

    var windowEl = $(window);
    var tp = window.innerHeight * 0.2 < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() ? __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 : window.innerHeight * 0.2;
    var header = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].header).render();
    var html = $.templates(__WEBPACK_IMPORTED_MODULE_0__view__["a" /* default */].body).render({ item: tmp });

    $(document.body).append(html);

    $('#ncit_details').dialog({
      modal: false,
      position: { my: "center top+" + tp, at: "center top", of: $('#docs-container') },
      width: 600,
      height: 450,
      minWidth: 420,
      maxWidth: 800,
      minHeight: 350,
      maxHeight: 650,
      open: function open() {
        //add new custom header
        $(this).prev('.ui-dialog-titlebar').css('padding-top', '3.5em').html(header);

        var target = $(this).parent();
        target.find('.ui-dialog-titlebar-close').html('');
        if (target.offset().top - windowEl.scrollTop() < __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset()) {
          target.css('top', windowEl.scrollTop() + __WEBPACK_IMPORTED_MODULE_2__shared__["a" /* default */].headerOffset() + 20 + 'px');
        }

        $('#close_ncit_details').bind('click', function () {
          $("#ncit_details").dialog('close');
        });
      },
      close: function close() {
        $(this).remove();
      }
    }).parent().draggable({
      containment: '#docs-container'
    });
  });
}

/***/ }),
/* 24 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
var tmpl = {
  header: '<div class="dialog__header">' + '<div class="dialog__titlebar">' + '<span id="ui-id-4" class="ui-dialog-title">NCIt Terms & Properties</span>' + '<button type="button" id="close_ncit_details" class="ui-button ui-corner-all ui-widget ui-button-icon-only ui-dialog-titlebar-close" title="Close"></button>' + '</div>' + '</div>',
  body: '<div id="ncit_details"><div class="ncit__content">' + '<p><b>Preferred Name:</b> {{:item.name}}</p>' + '{{if item.definition !== undefined}}' + '<p><b>Definition:</b> {{:item.definition}}</p>' + '{{/if}}' + '<p><b>NCI Thesaurus Code:</b>' + '<a href="https://ncit.nci.nih.gov/ncitbrowser/pages/concept_details.jsf?dictionary=NCI_Thesaurus&code={{:item.code}}" target="_blank"">{{:item.code}}</a>' + '</p>' + '{{if item.synonyms.length }}' + '<p><b>Synonyms &amp; Abbreviations:</b></p>' + '<p>{{for item.synonyms}}{{:}}</br>{{/for}}</p>' + '{{/if}}' + '<p><a href="https://ncit.nci.nih.gov/ncitbrowser/pages/concept_details.jsf?dictionary=NCI_Thesaurus&code={{:item.code }}" target="_blank"">more details</p>' + '</div>' + '</div>'
};

/* harmony default export */ __webpack_exports__["a"] = (tmpl);

/***/ }),
/* 25 */
/***/ (function(module, exports) {

// removed by extract-text-webpack-plugin

/***/ })
/******/ ]);