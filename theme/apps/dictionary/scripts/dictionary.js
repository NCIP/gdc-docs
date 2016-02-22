(function(_, d3, fetch, Promise) {
  'use strict';

  if (typeof window.Dictionary !== 'undefined') {
    throw new Error('Dictionary app looks to be already defined or is in use by someone else - Aborting....');
  }

  window.Dictionary = Dictionary;

  var _DICTIONARY_CONSTANTS = {
    VIEWS: {
      TABLE: 1
    },
    DICTIONARY_KEY_ORDER: [
      // clinical
      {'clinical': ['demographic', 'diagnosis', 'family_history', 'exposure', 'treatment']},

      // biospecimen
      {'biospecimen': ['sample', 'portion', 'analyte', 'aliquot']},

      // data bundle
      {'data_bundle': ['read_group', 'submitted_file',
        'clinical_data_bundle', 'biospecimen_data_bundle',
        'slide_data_bundle', 'pathology_data_bundle']},
      // annotation,
      {'annotation': ['annotation']}
    ],
    RESPONSE_TYPE: {
      TEXT: 'text/plain',
      JSON: 'application/json'
    },
    WEB_SERVICE: {
      DEFAULT_URL: '',
      CONTEXT_PATTERN: '/auth/api/v0/submission/${program}/${project}/_dictionary/${dictionary_name}',
      DEFAULT_PROGRAM: 'CGCI',
      DEFAULT_PROJECT: 'BLGSP',
      DEFAULT_DICTIONARY: '_all'
    }

  };

  var _defaultOptions = {
    onInitFn: _.noop,
    beforeRenderFn: _.noop,
    afterRenderFn: _.noop,
    dictionaryData: null, // if not null will use this as the data source
    dataSourceBaseHost: _DICTIONARY_CONSTANTS.WEB_SERVICE.DEFAULT_URL, // override internal host defaults
    dataSourceContextPattern: _DICTIONARY_CONSTANTS.WEB_SERVICE.CONTEXT_PATTERN,
    defaultView: _DICTIONARY_CONSTANTS.VIEWS.TABLE
  };


  ///////// Helpers
  function _checkResponseStatus(response) {
    if (response.status >= 200 && response.status < 300) {
      return response;
    }
    else {
      var error = new Error(response.statusText);
      error.response = response;
      throw error;
    }
  }

  function _responseParseText(response) {
    return response.text();
  }

  function _responseParseJSON(response) {
    return response.json();
  }

  function parseContextPattern(contextPattern, patternMapping) {

    var patternMapping = patternMapping || {},
        newContextPattern = contextPattern,
        patterns = ['program' , 'project', 'dictionary_name'];

    // Set pattern defaults if none is supplied...
    patternMapping.program = patternMapping.program || _DICTIONARY_CONSTANTS.WEB_SERVICE.DEFAULT_PROGRAM;
    patternMapping.project = patternMapping.project || _DICTIONARY_CONSTANTS.WEB_SERVICE.DEFAULT_PROJECT;
    patternMapping.dictionary_name =   patternMapping.dictionary_name || _DICTIONARY_CONSTANTS.WEB_SERVICE.DEFAULT_DICTIONARY;

    for (var i = 0; i < patterns.length; i++) {
      var pattern = patterns[i];

      // Warning matching is case sensitive!
      var tokenPattern = '${' + pattern + '}';

      if (contextPattern.indexOf(tokenPattern) >=0 && _.isString(patternMapping[pattern])) {
        // Replace all occurrences - commented out for now
        //newContextPattern = newContextPattern.split(tokenPattern).join(patternMapping[pattern]);
        newContextPattern = newContextPattern.replace(tokenPattern, patternMapping[pattern])
      }

    }

    return newContextPattern;

  }
  ////////


  function _initDictionaryData(data) {

    if (! data || _.isEmpty(data) || ! _.has(data, '_definitions')) {
      console.warn('Dictionary endpoint has returned invalid dictionary data! Dictionary data population aborting...', 'Data: ', data);
      return null;
    }

    var defData = data._definitions;

    delete data._definitions;

    var dictDataList = data,
        dictionaryData = {definitions:  defData, dictionaries: [], dictionaryMap: dictDataList, mapByCategory: {}};

    for (var dictionaryTitle in dictDataList) {
      var dictionary = dictDataList[dictionaryTitle],
          dictionaryCategory = dictionary.category || 'Unknown';

      if (dictDataList.hasOwnProperty(dictionaryTitle)) {
        dictionaryData.dictionaries.push(dictionary);
      }

      if (! _.isArray(dictionaryData.mapByCategory[dictionaryCategory]) ) {
        dictionaryData.mapByCategory[dictionaryCategory] = [];
      }

      dictionaryData.mapByCategory[dictionaryCategory].push(dictionary);
    }

    return dictionaryData;
  }

  /////////////////////////////////////////////////////////

  function Dictionary(targetEl, options) {
    var _dictionary = this;

    _dictionary._containerEl = null;
    _dictionary._options = _defaultOptions;
    _dictionary._data = null;


    if (_.isObject(options)) {
      _.assign(_dictionary._options, options);
    }

    _dictionary.containerElement = function(el) {
      // Check if the element is a dom node
      if ( Object.prototype.toString.call(el).toUpperCase().indexOf('HTML') >= 0 ) {
        _dictionary._containerEl = el;
      }

      return _dictionary._containerEl;
    };

    _dictionary.containerElement(targetEl || document.body);
    _dictionary.CONSTANTS = _DICTIONARY_CONSTANTS;

    _dictionary.init();
  }


  Dictionary.prototype.init = function(shouldRefresh) {

    var _dictionary = this;

    if (_dictionary._data === null || shouldRefresh === true) {
      _dictionary.getDataFromSource().then(function(rawDictionaryData) {
        _dictionary._data = _initDictionaryData(rawDictionaryData);
      });
    }

  };


  Dictionary.prototype.getDataFromSource = function(responseType) {
    var webServiceURL = this._options.dataSourceBaseHost + parseContextPattern(this._options.dataSourceContextPattern),
        fetchOptions = {
          //credentials: 'include',
          method: 'get'
        },
        responseParseFn = _responseParseJSON,
        responseMimeType = _.isString(responseType) ? responseType : '';

    switch(responseMimeType) {
      case _DICTIONARY_CONSTANTS.RESPONSE_TYPE.TEXT:
        responseParseFn = _responseParseText;
        break;
      default:
        break;
    }

    return new Promise(function(resolve, reject) {
      fetch(webServiceURL, fetchOptions)
        .then(_checkResponseStatus)
        .then(responseParseFn)
        .then(function (responseData) {
          //resolve original promise
          console.log('Request Succeeded - Data: ', responseData);
          resolve(responseData);
        })
        .catch(function (error) {
          console.log('Request Failed - Error: ', error);
          reject(error);
        });
    });
  };

  Dictionary.prototype.getData = function() {
    return this._data;
  };

  // Manual call to render self
  Dictionary.prototype.render = function() {

  };

  // Clean up
  Dictionary.prototype.destroy = function() {

  };



})(_.noConflict(), d3, window.fetch, window.Promise);