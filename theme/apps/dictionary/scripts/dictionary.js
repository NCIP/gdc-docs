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
    _dictionary._currentURLParams = _getParamsFromURL();

    _dictionary.containerElement(targetEl || document.body);
    _dictionary.CONSTANTS = _DICTIONARY_CONSTANTS;

    _dictionary.init();
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

    var _dictionary = this,
        _urlParams = _dictionary._currentURLParams;

    if (! _.isEmpty(_urlParams) ) {

      if (_.isString(_urlParams.view)) {
        _dictionary._currentView = _urlParams.view;
      }

    }

    if (_dictionary._data === null || shouldRefresh === true) {

      _fetchTemplate(_DICTIONARY_CONSTANTS.TEMPLATES.MAIN_DICTIONARY)
        .then(function (html) {
          _dictionary._containerEl.innerHTML = html;

          // Place D3 entry point
          _dictionary._d3Containers.rootEl = d3.select('#dictionary-inner-container');
        })
        .then(function () {
          _dictionary.getDataFromSource()
            .then(function (rawDictionaryData) {
              _dictionary._data = _initDictionaryData(rawDictionaryData);
              _dictionary._d3Containers.views = _getD3ViewsForDictionary(_dictionary._data, _dictionary.viewListener);
              _dictionary.setView(_dictionary._options.defaultView);
              _dictionary.render();

              _updatePageScroll(_dictionary._currentURLParams.anchor);
          });
        });
    }
  };

  Dictionary.prototype.viewListener = function(viewUpdateObj) {
    console.log(viewUpdateObj);
    console.log('view Listener invoked: ', viewUpdateObj.view.getState());

    switch( viewUpdateObj.eventType ) {

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INNER_NAV:

        if (_.isString(viewUpdateObj.params.id)) {
          _updatePageScroll(viewUpdateObj.params.id);
        }
        break;

      default:
        break;
    }
  };

  Dictionary.prototype.getDataFromSource = function(responseType) {
    var webServiceURL = this._options.dataSourceBaseHost + _parseContextPattern(this._options.dataSourceContextPattern);

    return _fetch(webServiceURL, responseType);
  };

  Dictionary.prototype.getSourceData = function() {
    return this._data;
  };

  Dictionary.prototype.setView = function(view) {
    var _dictionary = this;

    _dictionary._currentView = _dictionary._d3Containers.views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][view].view;

    return _dictionary;
  };

  // Manual call to render self
  Dictionary.prototype.render = function() {
    var _dictionary = this;

    _dictionary._currentView.show().render();


    return _dictionary;
  };

  // Clean up
  Dictionary.prototype.destroy = function() {
    console.log('Cleaning up the dictionary...');
  };

  /////////////////////////////////////////////////////////
  // Dictionary Views
  /////////////////////////////////////////////////////////

  function View(d3ContainerSelection, dictionaryData, actionCallbackFn) {
    var _view = this;

    _view._d3ContainerSelection = d3ContainerSelection;
    _view._dictionaryData = dictionaryData;
    _view._name = 'Unknown View';
    _view._callbackFn = actionCallbackFn || _.noop;
    _view._isHidden = true;

  }

  View.prototype.render = function() {
    var _view = this;

    console.log('Rendering!');

    _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.RENDERED;
    _view._callbackFn.call(null, new ViewUpdateObject(this));

    return _view;
  };

  View.prototype.show = function() {
    var _view = this;

    console.log('Showing!');

    _view._isHidden = false;
    _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.ENTER;
    _view._d3ContainerSelection.transition().duration(10).style('opacity', 1);

    _view._callbackFn.call(null, new ViewUpdateObject(this));

    return _view;
  };

  View.prototype.hide = function() {
    var _view = this;
    console.log('Hiding!');

    _view._isHidden = true;
    _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.EXIT;
    _view._d3ContainerSelection.transition().duration(10).style('opacity', 0);
    _view._callbackFn.call(null, new ViewUpdateObject(this));

    return _view;
  };

  View.prototype.getState = function() {
    return this._state;
  };

  /////////////////////////////////////////////////////////
  // TableEntityListView
  /////////////////////////////////////////////////////////
  var TableEntityListView = (function() {


    function TableEntityListView() {

      var _tableEntityListView = this;
      // Inherit from View
      View.apply(_tableEntityListView, arguments);


      _tableEntityListView._name = _DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST;


    }

    TableEntityListView.prototype = View.prototype;

    TableEntityListView.prototype.renderEntity = function(category, categoryData) {
      var _tableEntityListView = this;

      var entityTable = _tableEntityListView._d3ContainerSelection
        .append('table')
        .classed('dictionary-entity-table card', true)
        .attr('id', 'dictionary-entity-' + category);

      var tHead = entityTable.append('thead'),
          tBody = entityTable.append('tbody');

      tHead.append('tr')
        .append('th')
        .classed('dictionary-entity-header', true)
        .append('a')
        .attr('id', category)
        .attr('href',  '#?view=' + _tableEntityListView._name + '&anchor=' + category)
        .on('click', function() {
          _tableEntityListView._callbackFn.call(
            null, new ViewUpdateObject(_tableEntityListView,  _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INNER_NAV, {id: category})
          );
        })
        .html(function() {
          return '<i class="fa fa-book"></i> ' + _.get(_DICTIONARY_CONSTANTS.DICTIONARY_ENTITY_MAP, category.toLowerCase(), category);
        });

      var tRows = tBody.selectAll('tr')
            .data(categoryData)
            .enter()
            .append('tr');


      tRows.selectAll('td')
        .data(function(row) {
          return [{id: row.id,  title: row.title}];
        })
        .enter()
        .append('td')
        .classed('dictionary-entity-list-item', true)
        .append('a')
        .attr('title', function(data) {
          return 'View the definition of ' + data.title;
        })
        .attr('id', function(data) { return data.id; })
        .attr('href', function(data) {
          return '#?view=' + _tableEntityListView._name + '&id=' + data.id;
        })
        .on('click', function(data) {
          _tableEntityListView._callbackFn.call(
            null, new ViewUpdateObject(_tableEntityListView,  _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.NAV, {id: data.id})
          );
        })
        .text(function(data) { return data.title; });

      return entityTable;
    };

    TableEntityListView.prototype.render = function () {
      var _tableEntityListView = this;

     console.log(_tableEntityListView._dictionaryData);

      var categoryMap = _tableEntityListView._dictionaryData.dictionaryMapByCategory;

      for (var category in categoryMap) {
        if (categoryMap.hasOwnProperty(category)) {
          _tableEntityListView.renderEntity(category, categoryMap[category]);
        }
      }

      console.log('TableEntityListView Rendering!');





      _tableEntityListView._state = _DICTIONARY_CONSTANTS.VIEW_STATE.RENDERED;
      _tableEntityListView._callbackFn.call(null, new ViewUpdateObject(_tableEntityListView));
    };
    return TableEntityListView;
  })();

  /////////////////////////////////////////////////////////
  function TableDefinitionsView() {

    var _tableDefView = this;

    // Inherit from View
    View.apply(_tableDefView, arguments);


    _tableDefView._name = _DICTIONARY_CONSTANTS.VIEWS.TABLE.DEFINITION;

  }

  TableDefinitionsView.prototype = View.prototype;





  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////




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
        _ID: 'table',
        ENTITY_LIST: 'table-entity-list',
        DEFINITION: 'table-definition'
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
      tbd: 'To Be Determined...'
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
      DEFAULT_URL: 'https://gdc-api.nci.nih.gov',
      CONTEXT_PATTERN: '/auth/api/v0/submission/${program}/${project}/_dictionary/${dictionary_name}',
      DEFAULT_PROGRAM: 'CGCI',
      DEFAULT_PROJECT: 'BLGSP',
      DEFAULT_DICTIONARY: '_all'
    },
    TEMPLATES: {
      RELATIVE_DIR: '/html-shards',
      MAIN_DICTIONARY: 'dictionary.html'
    },
    BROWSER_CAPABILITIES: {
      SMOOTH_SCROLL: 'scrollBehavior' in document.documentElement.style
    }

  };

  ///////////////////////////////////////////////
  // Dictionary Default Options
  ///////////////////////////////////////////////
  var _defaultOptions = {
    onInitFn: _.noop,
    beforeRenderFn: _.noop,
    afterRenderFn: _.noop,
    dictionaryData: null, // if not null will use this as the data source
    dataSourceBaseHost: _DICTIONARY_CONSTANTS.WEB_SERVICE.DEFAULT_URL, // override internal host defaults
    dataSourceContextPattern: _DICTIONARY_CONSTANTS.WEB_SERVICE.CONTEXT_PATTERN,
    defaultView: _DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST
  };

  function ViewUpdateObject(viewObject, viewEventType, viewParams) {
    var _viewObject = this;

    _viewObject.view = viewObject;
    _viewObject.eventType = viewEventType || _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.DEFAULT;
    _viewObject.params = viewParams;
  }

  ///////////////////////////////////////////////
  // Initialize D3 Views
  ///////////////////////////////////////////////
  function _getD3ViewsForDictionary(dictionaryData, actionCallbackFn) {
    var views =  {},
        tableViews = {
          summary: d3.select('#dictionary-view-table-summary'),
          detailed: d3.select('#dictionary-view-table-detail')
        };

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID] = {};

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][_DICTIONARY_CONSTANTS.VIEWS.TABLE.ENTITY_LIST] = {
      el: tableViews.summary,
      view: new TableEntityListView(tableViews.summary, dictionaryData, actionCallbackFn)
    };

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][_DICTIONARY_CONSTANTS.VIEWS.TABLE.DEFINITION] = {
      el: tableViews.detailed,
      view: new TableDefinitionsView(tableViews.detailed, dictionaryData, actionCallbackFn)
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

  function _getParamsFromURL() {
    var view = null,
        hash = window.location.hash;

    if (! hash) {
      return view;
    }

    var argStr = hash.substr(2); // eat the hash and ? char

    var argTokens = argStr.split('&'),
        argParams = {};

    for (var i = 0; i < argTokens.length; i++) {
      var keyVals = argTokens[i].split('=');

      if (keyVals.length !== 2) {
        console.log('Unexpected hash arguments. Got: ', keyVals);
        continue;
      }

      argParams[keyVals[0]] = keyVals[1];
    }

    console.log(argParams);


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
      console.log('View found in ULR = ', view);
      argParams.view = view;
    }
    else {
      delete argParams.view;
    }

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

  function _fetch(url, responseType) {
    var webServiceURL = url,
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
  }

  function _fetchTemplate(templateFile) {
    return _fetch(  _DICTIONARY_CONSTANTS.APP_ABSOLUTE_DIR +
                    _DICTIONARY_CONSTANTS.TEMPLATES.RELATIVE_DIR + '/' +
                    templateFile, _DICTIONARY_CONSTANTS.RESPONSE_TYPE.TEXT  );
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

    var defData = data._definitions;

    delete data._definitions;

    var dictDataList = data,
      dictionaryData = {definitions:  defData, dictionaries: [], dictionaryMap: dictDataList, dictionaryMapByCategory: {}};

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