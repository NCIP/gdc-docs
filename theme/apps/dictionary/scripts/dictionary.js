(function(_, d3, fetch, Promise) {
  'use strict';

  /////////////////////////////////////////////////////////
  // Dictionary Main App
  /////////////////////////////////////////////////////////
  function Dictionary(targetEl, options) {
    var _dictionary = this;

    _dictionary._containerEl = null;
    _dictionary._d3Containers = {};
    _dictionary._options = _defaultOptions;
    _dictionary._data = null;



    if (_.isObject(options)) {
      _.assign(_dictionary._options, options);
    }

    _dictionary._currentViewMode = _dictionary._options.defaultView;
    _dictionary._currentView = null;

    _dictionary.containerElement(targetEl || document.body);
    _dictionary.CONSTANTS = _DICTIONARY_CONSTANTS;

    _dictionary.init();

    if (_DICTIONARY_CONSTANTS.BROWSER_CAPABILITIES.HASH_CHANGE_EVENT) {

      window.onhashchange = function (event) {

        if (window.location.href.toLowerCase().indexOf('dictionary/viewer') >= 0) {
          _dictionary.triggerViewEventFromURL(event);
        }

      };
    }
  }

  /////////////////////////////////////////////////////////
  // Get/Set the dictionary's parent container
  /////////////////////////////////////////////////////////
  Dictionary.prototype.containerElement  = function(el) {
    // Check if the element is a DOM node
    if ( Object.prototype.toString.call(el).toUpperCase().indexOf('HTML') >= 0 ) {
      this._containerEl = el;
    }

    return this._containerEl;
  };

  /////////////////////////////////////////////////////////
  // Dictionary Initializer
  /////////////////////////////////////////////////////////
  Dictionary.prototype.init = function(shouldRefresh) {

    if (typeof Dictionary._Views === 'undefined') {
      console.warn('Could not find the Dictionary Views - are you sure you included the views dictionary JS File.\nAborting Dictionary Init.');
      return;
    }

    var _dictionary = this;

    if (shouldRefresh === true) {
      _dictionary.destroy();
    }

    var urlParams = _getParamsFromURL();

    if (_.has(urlParams, 'view')) {
      _dictionary._currentViewMode = urlParams.view;
    }
    else {
      _dictionary._currentViewMode = _dictionary._options.defaultView;
    }


    if (_dictionary._data === null) {

      _fetchTemplate(_DICTIONARY_CONSTANTS.TEMPLATES.MAIN_DICTIONARY)
        .then(function (html) {
          _dictionary._containerEl.innerHTML = html;

          // Place D3 entry point
          _dictionary._d3Containers.rootSelection = d3.select('#dictionary-inner-container');

          _dictionary._d3Containers.breadcrumbs = {};
          _dictionary._d3Containers.breadcrumbs.breadcrumbSelection = d3.select('#dictionary-nav-breadcrumb');
          _dictionary._d3Containers.breadcrumbs.breadcrumbStackSelection = d3.select('#dictionary-nav-view-stack');

        })
        .then(function () {
          _dictionary.getDataFromSource()
            .then(function (rawDictionaryData) {

              _dictionary._data = _initDictionaryData(rawDictionaryData);

              _dictionary._d3Containers.views = _getD3ViewsForDictionary(_dictionary._data, _.bind(_dictionary.viewListener, _dictionary));

              _dictionary.triggerViewEventFromURL();
          });
        });
    }
  };

  Dictionary.prototype.viewListener = function(viewUpdateObj) {
    var _dictionary = this,
        view =  viewUpdateObj.view,
        params = viewUpdateObj.params;


    if (! view) {
      return;
    }

    console.log(viewUpdateObj);

    switch( viewUpdateObj.eventType ) {

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INNER_NAV:

        if (_.isString(params.id)) {
          _updatePageScroll(params.id);
        }
        break;

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.NAV:

        _dictionary._options.beforeRenderFn.call(_dictionary, _dictionary);

        var data = _dictionary._data,
            urlParams = _getParamsFromURL(),
            currentView = _dictionary.getCurrentView();

        if (currentView && currentView !== view) {
          currentView.hide();
        }

        if (_.has(params, 'id')) {
          data = {};
          _.assign(data, _dictionary._data.dictionaryMap[params.id], {terms: _.get(_dictionary._data, 'terms',{})});
        }

        _dictionary
          .setCurrentView(view.getViewName())
          .getCurrentView()
          .setDictionaryData(data)
          .show()
          .render();

        if (_.has(urlParams, 'anchor')) {
          _updatePageScroll(urlParams.anchor);
        }

        _dictionary.updateBreadcrumb();

        _dictionary._options.afterRenderFn.call(_dictionary, _dictionary);

        break;

      default:
        break;
    }
  };

  Dictionary.prototype.getCurrentViewName = function() {
    return this._currentViewMode;
  };

  Dictionary.prototype.getDataFromSource = function(responseType) {
    var webServiceURL = this._options.dataSourceBaseHost + _parseContextPattern(this._options.dataSourceContextPattern);

    return _fetch(webServiceURL, responseType);
  };


  Dictionary.prototype.getDictionaryTemplates = function(category, excludes, dataFormat) {
    var _dictionary = this,
        entityCategory = category || '',
        fileFormat = dataFormat || _DICTIONARY_CONSTANTS.END_POINT.ENDPOINT_PARAMS.TEMPLATE.TSV_TYPE,
        entityExclusions = _.isArray(excludes) && excludes.length ? excludes : [],
        params = {format: fileFormat},
        webServiceURL = this._options.dataSourceBaseHost + _parseContextPattern(_DICTIONARY_CONSTANTS.END_POINT.CONTEXT_TEMPLATE_PATTERN, {dictionary_name: ''}),
        containerEl = _dictionary._containerEl;

    if (entityExclusions.length) {
      params.exclude = entityExclusions.join(',');
    }

    if (entityCategory) {
      params.category = category;
    }

    var f = _createHiddenForm(containerEl, webServiceURL, params);
    f.submit();
    containerEl.removeChild(f);

  };

  Dictionary.prototype.getSourceData = function() {
    return this._data;
  };

  Dictionary.prototype._getView = function(viewName) {
    var _dictionary = this;

    if (! _dictionary._d3Containers.views ||
        ! _.has(_dictionary._d3Containers.views, _DICTIONARY_CONSTANTS.VIEWS.TABLE._ID) ||
        ! _.has(_dictionary._d3Containers.views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID], viewName)) {
      return null;
    }

    return _dictionary._d3Containers.views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][viewName].view;
  };

  Dictionary.prototype.setCurrentView = function(viewName) {
    var _dictionary = this;

    _dictionary._currentView = _dictionary._getView(viewName);
    _dictionary._currentViewMode = viewName;

    return _dictionary;
  };

  Dictionary.prototype.getCurrentView = function() {
    var _dictionary = this;

    return  _dictionary._currentView;
  };

  Dictionary.prototype.updateBreadcrumb = function() {
    var _dictionary = this,
        currentView = _dictionary.getCurrentView(),
        currentViewName = currentView.getViewName(),
        breadcrumbSelection = _dictionary._d3Containers.breadcrumbs.breadcrumbSelection,
        breadcrumbStack = _dictionary._d3Containers.breadcrumbs.breadcrumbStackSelection,
        breadcrumbName = currentView.getBreadcrumbName() || 'Unknown',
        styleOpts = {display: 'none'};


    if (currentViewName !== _DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST && breadcrumbName) {
      styleOpts.display = 'block';
    }

    if (breadcrumbSelection) {
      breadcrumbSelection.style(styleOpts);
    }


    breadcrumbStack.html('');

    breadcrumbStack.append('i').classed('fa fa-chevron-right', true);

    breadcrumbStack.append('span').text(' ' + breadcrumbName);

  };

  Dictionary.prototype.triggerViewEventFromURL = function(event) {
    var _dictionary = this,
        urlParams = _getParamsFromURL(true),
        viewMode = _dictionary._options.defaultView,
        view = null,
        params = null;

    if (event) {
      console.log(event);
      event.preventDefault();
      event.stopPropagation();
    }

    if (_.has(urlParams, 'view')) {
      viewMode = urlParams.view;
    }

    view = _dictionary._getView(viewMode);

    if (_.has(urlParams, 'id')) {
      params = {id: urlParams.id};
    }

    var viewEvent = new Dictionary._ViewUpdateObject(view, _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.NAV, params);

    _dictionary.viewListener(viewEvent);

  };


  // Manual call to render self
  Dictionary.prototype.render = function(renderParams) {
    var _dictionary = this;

    _dictionary._currentView.show().render(renderParams);


    return _dictionary;
  };

  // Clean up
  Dictionary.prototype.destroy = function() {
    console.log('Cleaning up the dictionary...');

    if (_DICTIONARY_CONSTANTS.BROWSER_CAPABILITIES.HASH_CHANGE_EVENT) {
      window.onhashchange = _.noop;
    }

    var _dictionary = this;

    _urlParamsCache = null;
    _dictionary._data = null;
  };


  /////////////////////////////////////////////////////////
  // Create a global access point to this Dictionary
  /////////////////////////////////////////////////////////
  if (typeof window.Dictionary !== 'undefined') {
    throw new Error('Dictionary app looks to be already defined or is in use by someone else - Aborting....');
  }

  window.Dictionary = Dictionary;

  ///////////////////////////////////////////////
  // Dictionary Constants
  ///////////////////////////////////////////////
  var _DICTIONARY_CONSTANTS = {
    APP_ABSOLUTE_DIR: '/apps/dictionary',
    VIEWS: {
      TABLE: {
        _ID: 'TABLE',
        ENTITY_LIST: 'table-entity-list',
        TERM_DEFINITION: 'table-definition-view'
      }
    },
    VIEW_STATE: {
      ENTER: 'enter', EXIT: 'exit', RENDERED: 'rendered'
    },
    VIEW_UPDATE_EVENT_TYPES: {
      DEFAULT:'update',
      NAV: 'nav',
      INNER_NAV: 'inner-nav'
    },
    DICTIONARY_ENTITY_MAP: {
      case: 'Case',
      clinical: 'Clinical',
      biospecimen: 'Biospecimen',
      data_bundle: 'Data Bundles',
      annotation: 'Annotation',
      data_file: 'Data Files',
      references: 'References',
      administrative: 'Administrative',
      tbd: 'References'
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
    END_POINT: {
      DEFAULT_URL: 'https://gdc-api.nci.nih.gov',
      //CONTEXT_PATTERN: '/auth/api/v0/submission/${program}/${project}/_dictionary/${dictionary_name}',
      CONTEXT_PROGRAM_PROJECT_PATTERN: '/v0/submission/${program}/${project}/_dictionary/${dictionary_name}',
      CONTEXT_PATTERN: '/v0/submission/_dictionary/${dictionary_name}',
      CONTEXT_TEMPLATE_PATTERN: '/v0/submission/template/${dictionary_name}',
      DEFAULT_PROGRAM: 'CGCI',
      DEFAULT_PROJECT: 'BLGSP',
      DEFAULT_DICTIONARY: '_all',
      ENDPOINT_PARAMS: {
        TEMPLATE: {
          JSON_TYPE: 'json',
          CSV_TYPE: 'csv',
          TSV_TYPE: 'tsv'
        }
      }
    },
    TEMPLATES: {
      RELATIVE_DIR: '/html-shards',
      MAIN_DICTIONARY: 'dictionary.html'
    },
    BROWSER_CAPABILITIES: {
      SMOOTH_SCROLL: 'scrollBehavior' in document.documentElement.style,
      HASH_CHANGE_EVENT: 'onhashchange' in window
    },
    DATA_FORMATS: {
      MISSING_VAL: '--',
      TEXT_FORMAT: 'text/plain',
      JSON_FORMAT: 'application/json',
      TSV_FORMAT: 'text/tab-separated-values',
      XML_FORMAT: 'text/xml',
      CSV_FORMAT: 'text/csv',
      BLOB_FORMAT: 'octet/stream'
    }

  };

  // Memoize constants so the views can see them
  Dictionary._DICTIONARY_CONSTANTS = _DICTIONARY_CONSTANTS;
  Dictionary._ = _;

  ///////////////////////////////////////////////
  // Dictionary Default Options
  ///////////////////////////////////////////////
  var _urlParamsCache = null;

  var _defaultOptions = {
    onInitFn: _.noop,
    beforeRenderFn: _.noop,
    afterRenderFn: _.noop,
    dictionaryData: null, // if not null will use this as the data source
    dataSourceBaseHost: _DICTIONARY_CONSTANTS.END_POINT.DEFAULT_URL, // override internal host defaults
    dataSourceContextPattern: _DICTIONARY_CONSTANTS.END_POINT.CONTEXT_PATTERN,
    defaultView: _DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST
  };


  ///////////////////////////////////////////////
  // Initialize D3 Views
  ///////////////////////////////////////////////
  function _getD3ViewsForDictionary(dictionaryData, actionCallbackFn) {
    var views =  {},
        tableViews = {
          summary: d3.select('#dictionary-view-table-summary'),
          detailed: d3.select('#dictionary-view-table-detail')
        },
        urlParams = _getParamsFromURL();

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID] = {};

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][_DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST] = {
      el: tableViews.summary,
      view: new  Dictionary._Views.TableEntityListView(tableViews.summary, dictionaryData, actionCallbackFn)
    };

    if (_.has(urlParams, 'id')) {
      dictionaryData = dictionaryData.dictionaryMap[urlParams.id] || null;
    }

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][_DICTIONARY_CONSTANTS.VIEWS.TABLE.TERM_DEFINITION] = {
      el: tableViews.detailed,
      view: new  Dictionary._Views.TableDefinitionsView(tableViews.detailed, dictionaryData, actionCallbackFn)
    };


    return views;
  }

  function _updatePageScroll(anchor) {
    if (_.isString(anchor)) {
      _scrollTo(anchor);
    }
  }

  function _scrollTo(id) {
    var isSmoothScrollSupported = _DICTIONARY_CONSTANTS.BROWSER_CAPABILITIES.SMOOTH_SCROLL,
        node = document.getElementById(id);

    if (node) {
      var offset = node.getBoundingClientRect();

      var options = {
        behavior: 'smooth',
        left: 0,
        // Calculate the absolute offset and compensate for the menu bar
        top: Math.max(0, offset.top - 80) + window.scrollY
      };

      if (isSmoothScrollSupported) {
        // Native smooth scrolling
        window.scrollTo(options);
      }
      else {
        // Old way scrolling without effects
        window.scrollTo(options.left, options.top);
      }
    }

  }

  function _createHiddenForm(parentEl, action, params, method) {
    var formMethod = method || 'GET';

    var _form_ = document.createElement('form');

    _form_.setAttribute('method', formMethod);
    _form_.setAttribute('action', action);

    _form_.style.display = 'none';

    if (_.isObject(params)) {

      var paramKeys = _.keys(params);

      for (var i = 0; i < paramKeys.length; i++) {

        var paramName = paramKeys[i],
            input = document.createElement('input');

        input.setAttribute('type','hidden');
        input.setAttribute('name', paramName);
        input.setAttribute('value', params[paramName]);
        _form_.appendChild(input);

      }
    }

    parentEl.appendChild(_form_);

    return _form_;
  }

  function _getParamsFromURL(shouldNotCacheParams) {

    if (_urlParamsCache && shouldNotCacheParams !== true) {
      return _urlParamsCache;
    }

    var view = null,
        hash = window.location.hash,
        argParams = {};

    if (! hash) {
      return view;
    }

    if (hash.length <= 2) {
      return argParams;
    }

    var argStr = hash.substr(2); // eat the hash and ? char

    var argTokens = argStr.split('&');


    for (var i = 0; i < argTokens.length; i++) {
      var keyVals = argTokens[i].split('=');

      if (keyVals.length !== 2) {
        console.log('Unexpected hash arguments. Got: ', keyVals);
        continue;
      }

      argParams[keyVals[0]] = keyVals[1];
    }

    // Parse Logic
    //#view=BLAH&id=ANCHOR_NAME
    if (argParams.view) {

      for (var parentView in _DICTIONARY_CONSTANTS.VIEWS) {

        if (!_DICTIONARY_CONSTANTS.VIEWS.hasOwnProperty(parentView)) {
          continue;
        }

        var subViews = _DICTIONARY_CONSTANTS.VIEWS[parentView];

        for (var v in subViews) {
          if (v[0] === '_') {
            continue;
          }

          if (subViews[v] === argParams['view']) {
            view = subViews[v];
          }

        }
      }
    }

    if (view) {
      console.log('View found in URL = ', view);
      argParams.view = view;
    }
    else {
      delete argParams.view;
    }

    _urlParamsCache = argParams;

    return argParams;
  }

  ///////////////////////////////////////////////
  // Request/Response Helper Methods
  ///////////////////////////////////////////////
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

  function _parseContextPattern(contextPattern, patternMapping) {

    var patternMapping = patternMapping || {},
      newContextPattern = contextPattern,
      patterns = ['program' , 'project', 'dictionary_name'];

    // Set pattern defaults if none is supplied...
    patternMapping.program = patternMapping.program || _DICTIONARY_CONSTANTS.END_POINT.DEFAULT_PROGRAM;
    patternMapping.project = patternMapping.project || _DICTIONARY_CONSTANTS.END_POINT.DEFAULT_PROJECT;
    patternMapping.dictionary_name = typeof patternMapping.dictionary_name === 'string' ? patternMapping.dictionary_name : _DICTIONARY_CONSTANTS.END_POINT.DEFAULT_DICTIONARY;

    for (var i = 0; i < patterns.length; i++) {
      var pattern = patterns[i];

      // Warning matching is case sensitive!
      var tokenPattern = '${' + pattern + '}';

      if (contextPattern.indexOf(tokenPattern) >=0 && _.isString(patternMapping[pattern])) {
        // Replace all occurrences - commented out for now
        // ContextPattern = newContextPattern.split(tokenPattern).join(patternMapping[pattern]);
        newContextPattern = newContextPattern.replace(tokenPattern, patternMapping[pattern])
      }

    }

    return newContextPattern;

  }

  function _fetch(url, responseType) {
    var webServiceURL = url,
      fetchOptions = {
        //credentials: 'include',
        method: 'get'
      },
      responseParseFn = _responseParseJSON,
      responseMimeType = _.isString(responseType) ? responseType : '';

    switch(responseMimeType) {
      case _DICTIONARY_CONSTANTS.DATA_FORMATS.TEXT_FORMAT:
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
          console.warn('Request Succeeded - Data: ', responseData);
          resolve(responseData);
        })
        .catch(function (error) {
          console.warn('Request Failed - Error: ', error);
          reject(error);
        });
    });
  }

  function _fetchTemplate(templateFile) {
    return _fetch(  _DICTIONARY_CONSTANTS.APP_ABSOLUTE_DIR +
                    _DICTIONARY_CONSTANTS.TEMPLATES.RELATIVE_DIR + '/' +
                    templateFile, _DICTIONARY_CONSTANTS.DATA_FORMATS.TEXT_FORMAT  );
  }

  ///////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // Put the dictionary data returned into a useful format for later...
  ////////////////////////////////////////////////////////////////////////
  function _initDictionaryData(data) {

    if (! data || _.isEmpty(data) || ! _.has(data, '_definitions')) {
      console.warn('Dictionary endpoint has returned invalid dictionary data! Dictionary data population aborting...',
        'Received Data: ', data);
      return null;
    }

    var dictDataList = data,
        dictionaryData = {dictionaries: [], dictionaryMap: dictDataList, dictionaryMapByCategory: {}};

    var dictionaryKeys = _.keys(dictDataList);

    for (var i = 0; i < dictionaryKeys.length; i++) {
      var dictionaryKey = dictionaryKeys[i];

      // Add special private data prefixed with '_' to the dictionary data object directly...
      if (dictionaryKey[0] === '_') {
        dictionaryData[dictionaryKey.substr(1)] = dictDataList[dictionaryKey];
        delete dictDataList[dictionaryKey];
      }

    }

    console.log(dictionaryData);

    // Build our data structures and corresponding caches
    for (var dictionaryTitle in dictDataList) {
      var dictionary = dictDataList[dictionaryTitle],
        dictionaryCategory = dictionary.category || 'Unknown';

      if (dictDataList.hasOwnProperty(dictionaryTitle)) {
        dictionaryData.dictionaries.push(dictionary);
      }

      if (! _.isArray(dictionaryData.dictionaryMapByCategory[dictionaryCategory]) ) {
        dictionaryData.dictionaryMapByCategory[dictionaryCategory] = [];
      }

      dictionaryData.dictionaryMapByCategory[dictionaryCategory].push(dictionary);
    }

    return dictionaryData;
  }
  /////////////////////////////////////////////////////////

})(_.noConflict(), d3, window.fetch, window.Promise);