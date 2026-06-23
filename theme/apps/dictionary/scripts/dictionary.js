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

              _dictionary._d3Containers.views = _getD3ViewsForDictionary(_dictionary, _.bind(_dictionary.viewListener, _dictionary));

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

    // replace '^2' with superscript
    $(".property-description:contains('^')").each(function () {
      var el = $(this);
      var content = el.text();
      var matched = content.match(/\^\d*/);
      var toReplace = content.replace(
        /\^\d*/,
        '<sup>' + matched[0].split('').slice(1, Infinity).join('') + '</sup>'
      );
      el.html(toReplace);
    });

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

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.TEMPLATE_DOWNLOAD_BY_CATEGORY_REQUESTED:

        if (_.has(params, 'id')) {
         _dictionary.getDictionaryTemplates(params.id, _.get(params, 'excludes', null));
        }

        break;

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.TEMPLATE_DOWNLOAD_REQUESTED:

        if (_.has(params, 'id')) {
          _dictionary.getDictionaryTemplate(params.id);
        }

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

  Dictionary.prototype.getDictionaryTemplateURL = function(dictionaryID, dataFormat) {
    var _dictionary = this,
      fileFormat = dataFormat || _dictionary._options.defaultTemplateDownloadFormat,
      webServiceURL = _dictionary._options.dataSourceBaseHost + _parseContextPattern(_DICTIONARY_CONSTANTS.END_POINT.CONTEXT_TEMPLATE_PATTERN, {dictionary_name: dictionaryID});

    return webServiceURL + '?format=' + fileFormat;
  };

  Dictionary.prototype.getDictionaryTemplate = function(dictionaryID, dataFormat) {
    var _dictionary = this,
        fileFormat = dataFormat || _dictionary._options.defaultTemplateDownloadFormat,
        params = {format: fileFormat},
        webServiceURL = _dictionary._options.dataSourceBaseHost + _parseContextPattern(_DICTIONARY_CONSTANTS.END_POINT.CONTEXT_TEMPLATE_PATTERN, {dictionary_name: dictionaryID}),
        containerEl = _dictionary._containerEl;

    var f = _createHiddenForm(containerEl, webServiceURL, params);
    f.submit();
    containerEl.removeChild(f);
  };

  Dictionary.prototype.getDictionaryTemplates = function(categories, excludes, dataFormat) {
    var _dictionary = this,
        entityCategories = categories,
        fileFormat = dataFormat || _dictionary._options.defaultTemplateDownloadFormat,
        entityExclusions = _.isArray(excludes) && excludes.length ? excludes : [],
        params = { format: fileFormat },
        webServiceURL = _dictionary._options.dataSourceBaseHost + _parseContextPattern(_DICTIONARY_CONSTANTS.END_POINT.CONTEXT_TEMPLATE_PATTERN, {dictionary_name: ''}),
        containerEl = _dictionary._containerEl;

    if (entityExclusions.length > 0) {
      console.warn('Excluding Entities: ', entityExclusions);
      params.exclude = entityExclusions.join(',');
    }

    if (_.isString(entityCategories)) {
      entityCategories = [entityCategories];
    }

    if (_.isArray(entityCategories) && entityCategories.length > 0) {
      params.categories = entityCategories.join(',');
      console.warn('Including Categories: ', entityCategories);
    }

    var f = _createHiddenForm(containerEl, webServiceURL, params);
    f.submit();
    containerEl.removeChild(f);

  };

  Dictionary.prototype.setDefaultDictionaryTemplateDownloadFormat = function(downloadFormat) {
    this._options.defaultTemplateDownloadFormat = downloadFormat;
  };

  Dictionary.prototype.getDefaultDictionaryTemplateDownloadFormat = function() {
    return this._options.defaultTemplateDownloadFormat;
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

  Dictionary.prototype.fetchDictionaryTemplate = function(templateFile) {
    return _fetchTemplate(templateFile);
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
      params = { id: urlParams.id };
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
      },
      _STATIC: {
        DICTIONARY_CONTROLS: '#dictionary-view-summary-controls'
      }
    },
    VIEW_STATE: {
      ENTER: 'enter', EXIT: 'exit', RENDERED: 'rendered'
    },
    VIEW_UPDATE_EVENT_TYPES: {
      DEFAULT:'update',
      NAV: 'nav',
      INNER_NAV: 'inner-nav',
      TEMPLATE_DOWNLOAD_BY_CATEGORY_REQUESTED: 'template-category-requested',
      TEMPLATE_DOWNLOAD_REQUESTED: 'template-dictionary-requested'
    },
    DICTIONARY_ENTITY_MAP: {
      case: 'Case',
      clinical: 'Clinical',
      biospecimen: 'Biospecimen',
      annotation: 'Annotations',
      data_file: 'Data Files',
      generated_data_file: 'Generated Data Files',
      references: 'References',
      administrative: 'Administrative',
      metadata_file: 'Metadata Files',
      tbd: 'References',
      index_file: 'Index',
      submittable_data_file: 'Submittable Data Files',
    },
    ENTITY_LIST_DICTIONARY_KEY_ORDER: ['case', 'clinical', 'biospecimen', 'submittable_data_file', 'generated_data_file', 'annotation', 'administrative', 'analysis', 'notation'],
    CATEGORY_TEMPLATE_DOWNLOAD_BLACKLIST: ['tbd', 'administrative', 'index_file', 'analysis', 'notation', 'generated_data_file'],
    CATEGORY_EXCLUDES: ['TBD'],
    CATEGORY_TEMPLATE_EXCLUDES: {
      clinical: ['clinical', 'clinical_test'],
      annotation: ['analysis', 'archive', 'publication', 'slide']
    },
    LINK_EXCLUDES: ['file', 'archive', 'clinical_test', 'submitted_methylation_beta_value'],
    PROPERTY_EXCLUDES: ['type', 'clinical_data_bundles', 'biospecimen_data_bundles', 'pathology_data_bundles', 'batch_id'],
    CATEGORY_TEMPLATE_INCLUDES: {
    },
    END_POINT: {
      DEFAULT_URL: 'https://api.gdc.cancer.gov', // TODO: env variable? 'http://localhost:5000'
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
    defaultView: _DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST,
    defaultTemplateDownloadFormat: _DICTIONARY_CONSTANTS.END_POINT.ENDPOINT_PARAMS.TEMPLATE.TSV_TYPE
  };


  ///////////////////////////////////////////////
  // Initialize D3 Views
  ///////////////////////////////////////////////
  function _getD3ViewsForDictionary(dictionary, actionCallbackFn) {
    var views =  {},
        tableViews = {
          summary: d3.select('#dictionary-view-table-summary'),
          detailed: d3.select('#dictionary-view-table-detail')
        },
        urlParams = _getParamsFromURL(),
        dictionaryData = dictionary._data;

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID] = {};

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][_DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST] = {
      el: tableViews.summary,
      view: new  Dictionary._Views.TableEntityListView(tableViews.summary, dictionary, actionCallbackFn, dictionary)
    };

    if (_.has(urlParams, 'id')) {
      dictionaryData = dictionaryData.dictionaryMap[urlParams.id] || null;
    }

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][_DICTIONARY_CONSTANTS.VIEWS.TABLE.TERM_DEFINITION] = {
      el: tableViews.detailed,
      view: new Dictionary._Views.TableDefinitionsView(tableViews.detailed, dictionaryData, actionCallbackFn, dictionary)
    };

    return views;
  }

  function _updatePageScroll(anchor) {
    if (_.isString(anchor)) {
      setTimeout(function() { _scrollTo(anchor); }, 50);
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

  var _fetchTemplateCache = {};

  function _fetchTemplate(templateFile, forceReload) {
    var cachedTemplate = _.get(_fetchTemplateCache, templateFile, false);

    if (_.isString(cachedTemplate) && forceReload !== true) {
      return new Promise(function(resolve) { resolve(cachedTemplate); });
    }

    var promise = _fetch(  _DICTIONARY_CONSTANTS.APP_ABSOLUTE_DIR +
                    _DICTIONARY_CONSTANTS.TEMPLATES.RELATIVE_DIR + '/' +
                    templateFile, _DICTIONARY_CONSTANTS.DATA_FORMATS.TEXT_FORMAT  );

    promise.then(function(html) {
      _fetchTemplateCache[templateFile] = html;
    });

    return promise;
  }

  ///////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // Put the dictionary data returned into a useful format for later...
  ////////////////////////////////////////////////////////////////////////
  function _sanitizeDictionaryData(d) {
    var data = _.cloneDeep(d);
    var dictionaryData = {};

    // TODO: Cleanup below hardcoding but unfortunately now this is necessary
    delete data.clinical;

    var dictionaryKeys = _.keys(data);

    for (var i = 0; i < dictionaryKeys.length; i++) {
      var dictionaryName = dictionaryKeys[i];
      var dictionary = data[dictionaryName];
      var dictionaryCategory = _.get(dictionary, 'category');

      if (dictionaryName === 'case') {
        dictionary.category = 'case';
      }

      if (dictionaryName === 'annotation') {
        dictionary.category = 'annotation';
      }

      dictionaryData[dictionaryName] = dictionary;

    }

    return dictionaryData;
  }

  function _initDictionaryData(data) {

    if (! data || _.isEmpty(data) || ! _.has(data, '_definitions')) {
      console.warn('Dictionary endpoint has returned invalid dictionary data! Dictionary data population aborting...',
        'Received Data: ', data);
      return null;
    }

    var dictDataList = _sanitizeDictionaryData(data);

    // Build our data structures and corresponding caches
    var r = Object.keys(dictDataList).reduce(function (acc, dictionaryTitle) {
      var unfilteredDict = dictDataList[dictionaryTitle];

      // filter out deprecated enums and properties but not required and preferred
      var dictionary = Object.assign(
        unfilteredDict,
        unfilteredDict.properties ? {
          properties: Object.keys(unfilteredDict.properties).reduce(
            (acc, key) => Object.assign(
              acc,
              ((unfilteredDict.deprecated || []).includes(key)
                && !([].concat(unfilteredDict.required || [], unfilteredDict.preferred || [])).includes(key)) ?
              {} : {[key]: Object.assign(
                unfilteredDict.properties[key],
                unfilteredDict.properties[key].enum ?
                {
                  enum: unfilteredDict.properties[key].enum.filter(
                    en => !(unfilteredDict.properties[key].deprecated_enum || [])
                    .includes(en)
                  )
                } : {}
              )}),
            {})
        } : {}
      );

      if (dictionaryTitle[0] === '_') {
        // Add special private data prefixed with '_' to the dictionary data object top level
        acc[dictionaryTitle.substr(1)] = dictionary;
      } else if (dictionary.category) {
        // otherwise add it to dictionaries array
        if (dictDataList.hasOwnProperty(dictionaryTitle)) {
          acc.dictionaries = acc.dictionaries.concat(dictionary);
        }
        // categorize by 'category' field except for data_file and metadata_file
        // furthur break those down into 'Submittable Data Files' (submittable) and 'Generated Data Files' (!submittable)
        // and add a ui_category key to reflect this (affects download template button display)
        dictionary.ui_category = dictionary.category;
        if (dictionary.category === 'data_file' || dictionary.category === 'metadata_file') {
          if (dictionary.submittable) {
            dictionary.ui_category = 'sumbittable_data_file';
            acc.dictionaryMapByCategory.submittable_data_file = acc.dictionaryMapByCategory.submittable_data_file.concat(dictionary);
          } else {
            dictionary.ui_category = 'generated_data_file';
            acc.dictionaryMapByCategory.generated_data_file = acc.dictionaryMapByCategory.generated_data_file.concat(dictionary);
          }
        } else {
          acc.dictionaryMapByCategory[dictionary.category] = (acc.dictionaryMapByCategory[dictionary.category] || []).concat(dictionary);
        }
      }
      return acc;
    }, { dictionaries: [], dictionaryMap: dictDataList, dictionaryMapByCategory: {'submittable_data_file': [], 'generated_data_file': [] } });
    return r;
  }
  /////////////////////////////////////////////////////////

})(_.noConflict(), d3, window.fetch, window.Promise);
