(function(_, d3) {
  'use strict';

  if (typeof window.Dictionary !== 'undefined') {
    throw new Error('Dictionary app looks to be already defined or is in use by someone else - Aborting....');
  }

  window.Dictionary = Dictionary;

  var _defaultOptions = {
    containerEl: document.body,
    onInitFn: _.noop(),
    onHideFn: _.noop(),
    dictionaryData: null
  };

  var CONSTANTS = {
    VIEWS: {
      'table': 1
    }
  };

  /////////////////////////////////////////////////////////

  function Dictionary(options) {
    var _dictionary = this;

    _dictionary.options = _defaultOptions;

    if (options) {
      _.extend(this.options, options)
    }

    _dictionary.CONSTANTS = CONSTANTS;

  }


  Dictionary.prototype.init = function(shouldRefresh) {

  };

  // Manual call to render self
  Dictionary.prototype.render = function() {

  };

  // Clean up
  Dictionary.prototype.destroy = function() {

  };



})(_.noConflict(), d3);