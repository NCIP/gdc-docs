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

    window.onpopstate = function(event) {
      if (window.location.href.toLowerCase().indexOf('dictionary/viewer') >= 0) {
        event.preventDefault();
        _dictionary.createViewEventFromURL();
      }

    };
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

    console.log( _dictionary._currentViewMode);

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

              /*_dictionary
                .setCurrentView(_dictionary._currentViewMode)
                .render();

              _updatePageScroll(urlParams.anchor);
              _dictionary.updateBreadcrumb();*/
              _dictionary.createViewEventFromURL();
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
    console.log('view Listener invoked: ', view.getState());

    switch( viewUpdateObj.eventType ) {

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INNER_NAV:

        if (_.isString(params.id)) {
          _updatePageScroll(params.id);
        }
        break;

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.NAV:

        var subViews = _.get(_DICTIONARY_CONSTANTS.VIEWS, view.getParentViewName(), false);

        if (view.getParentViewName() && subViews) {
          view.hide();

          if (view.getViewName() === subViews.ENTITY_LIST && _.has(params, 'id')) {

            _dictionary
              .setCurrentView(subViews.TERM_DEFINITION)
              .getCurrentView()
              .setDictionaryData(_dictionary._data.dictionaryMap[params.id])
              .show()
              .render();
          }

          _dictionary.updateBreadcrumb();
        }

        break;

      case _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INIT:

        var data = _dictionary._data,
            urlParams = _getParamsFromURL();

        if (_.has(params, 'id')) {
          data = _dictionary._data.dictionaryMap[params.id];
        }

        _dictionary
          .setCurrentView(view.getViewName())
          .getCurrentView(data)
          .setDictionaryData(_dictionary._data)
          .show()
          .render();

        if (_.has(urlParams, 'anchor')) {
          _updatePageScroll(urlParams.anchor);
        }

        _dictionary.updateBreadcrumb();

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

  Dictionary.prototype._getView = function(viewName) {
    var _dictionary = this;

    return _dictionary._d3Containers.views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][viewName].view;
  };

  Dictionary.prototype.setCurrentView = function(viewName) {
    var _dictionary = this;

    _dictionary._currentView = _dictionary._getView(viewName);

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

  Dictionary.prototype.createViewEventFromURL = function() {
    var _dictionary = this,
        urlParams = _getParamsFromURL(),
        viewMode = _dictionary._options.defaultView,
        view = null,
        params = null;

    if (_.has(urlParams, 'view')) {
      viewMode = urlParams.view;
    }

    view = _dictionary._getView(viewMode);

    if (_.has(urlParams, 'id')) {
      params = {id: urlParams.id};
    }

    _dictionary.viewListener(new ViewUpdateObject(view, _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INIT, params));

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

    var _dictionary = this;

    _urlParamsCache = null;
    _dictionary._data = null;
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
    _view._parentViewName = '';
    _view._prettyName = _view._name;
    _view._breadcrumbName = null;

    /////////////////////////////////////////////////////////
    // Public View API
    /////////////////////////////////////////////////////////
    _view.render = function() {
      console.log('Rendering!');

      _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.RENDERED;
      _view._callbackFn.call(null, new ViewUpdateObject(this));

      return _view;
    };

    _view.show = function() {
      console.log('Showing!');

      _view._isHidden = false;
      _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.ENTER;
      _view._d3ContainerSelection.style('display', 'block').transition().duration(250).style('opacity', 1);

      _view._callbackFn.call(null, new ViewUpdateObject(this));

      return _view;
    };

    _view.hide = function() {
      console.log('Hiding!');

      _view._isHidden = true;
      _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.EXIT;
      _view._d3ContainerSelection.transition().duration(250).style('opacity', 0);
      setTimeout(function() { _view._d3ContainerSelection.style('display', 'none'); }, 250);


      _view._callbackFn.call(null, new ViewUpdateObject(this));

      return _view;
    };

    _view.getState = function() {
      return _view._state;
    };

    _view.getViewName = function() {
      return _view._name;
    };

    _view.getParentViewName = function() {
      return _view._parentViewName;
    };

    _view.getPrettyName = function() {
      return _view._prettyName;
    };

    _view.setDictionaryData = function(data) {
      _view._dictionaryData = data;
      return _view;
    }

    _view.getBreadcrumbName = function() {
      return _view._breadcrumbName;
    };


  }


  /////////////////////////////////////////////////////////
  // TableEntityListView
  /////////////////////////////////////////////////////////
  var TableEntityListView = (function() {


    function TableEntityListView() {

      var _tableEntityListView = this;
      // Inherit from View
      View.apply(_tableEntityListView, arguments);

      _tableEntityListView._parentViewName = _DICTIONARY_CONSTANTS.VIEWS.TABLE._ID;
      _tableEntityListView._name = _DICTIONARY_CONSTANTS.VIEWS[_tableEntityListView._parentViewName].ENTITY_LIST;
      _tableEntityListView._prettyName = 'Entity List';

      /////////////////////////////////////////
      // Custom Render
      _tableEntityListView.render = function () {

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

    }

    // Inherit our defaults for View
    TableEntityListView.prototype = new View();

    TableEntityListView.prototype.renderEntity = function(category, categoryData) {
      var _tableEntityListView = this;

      var entityTable = _tableEntityListView._d3ContainerSelection
        .append('table')
        .classed('dictionary-entity-table card', true)
        .attr('id', 'dictionary-entity-' + category);

      var tHead = entityTable.append('thead'),
          tBody = entityTable.append('tbody');

      var getTooltipText = function() {

        var tooltipText = null;

        switch(_.first(categoryData).category) {
          case 'clinical':
            tooltipText = 'Cases must be registered in GDC before clinical, biospecimen, experiment and annotation data can be submitted.';
            break;
          default:
            break;
        }

        return tooltipText;
      };

      tHead.append('tr')
        .append('th')
        .attr('colspan', 2)
        .classed('dictionary-entity-header', true)
        .append('a')
        .classed('dictionary-tooltip', function() {
            return _.isString(getTooltipText());
        })
        .attr('id', category)
        .attr('href',  '#?view=' + _tableEntityListView._name + '&anchor=' + category)
        .on('click', function() {
          _tableEntityListView._callbackFn.call(
            null, new ViewUpdateObject(_tableEntityListView,  _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INNER_NAV, {id: category})
          );
        })
        .html(function() {
          var tooltipText = getTooltipText();
          return '<i class="fa fa-book"></i> ' + _.get(_DICTIONARY_CONSTANTS.DICTIONARY_ENTITY_MAP, category.toLowerCase(), category) +
                 (_.isString(tooltipText) ? '<span><i></i>' + tooltipText + '</span>' : '');
        });

      var tRows = tBody.selectAll('tr')
            .data(categoryData)
            .enter()
            .append('tr')
            .classed('dictionary-entity-list-item', true);


      tRows.selectAll('td')
        .data(function(row) {
          return [
            {id: row.id,  title: row.title, description: row.description},
            {id: row.id,  title: row.title, description: row.description}
          ];
        })
        .enter()
        .append('td')
        .classed('link', function(data, i) {
          var isLink = false;

          if (i === 0) {
            isLink = true;
          }
          return isLink;
        })
        .append('a')
        .attr('title', function(data) {
          return 'View details about ' + data.title;
        })
        .attr('id', function(data) { return data.id; })
        .attr('href', function(data) {
          return '#?view=' + _DICTIONARY_CONSTANTS.VIEWS[_tableEntityListView.getParentViewName()].TERM_DEFINITION + '&id=' + data.id;
        })
        .on('click', function(data) {
          _tableEntityListView._callbackFn.call(
            null, new ViewUpdateObject(_tableEntityListView,  _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.NAV, {id: data.id})
          );
        })
        .text(function(data, i) {
          var item = data.title;

          if (i === 1) {
            item = data.description;
          }

          return item;
        });

      return entityTable;
    };


    return TableEntityListView;
  })();

  /////////////////////////////////////////////////////////
  // TableDefinitionsView
  /////////////////////////////////////////////////////////
  var TableDefinitionsView = (function() {
    function TableDefinitionsView() {

      var _tableDefinitionView = this;

      // Inherit from View
      View.apply(_tableDefinitionView, arguments);

      _tableDefinitionView._parentViewName = _DICTIONARY_CONSTANTS.VIEWS.TABLE._ID;
      _tableDefinitionView._name = _DICTIONARY_CONSTANTS.VIEWS[_tableDefinitionView._parentViewName].TERM_DEFINITION;
      _tableDefinitionView._prettyName = 'Definition View';

      ////////////////////////////////////////////////////////////
      // Custom Render
      ////////////////////////////////////////////////////////////
      _tableDefinitionView.render = function() {

       if (_tableDefinitionView._dictionaryData) {
         _tableDefinitionView.renderDefinitionView(_tableDefinitionView._dictionaryData);
       }

        console.log('TableDefinitionsView Rendering!');

        _tableDefinitionView._state = _DICTIONARY_CONSTANTS.VIEW_STATE.RENDERED;
        _tableDefinitionView._callbackFn.call(null, new ViewUpdateObject(_tableDefinitionView));
      };

      _tableDefinitionView.setDictionaryData = function(data) {

        _tableDefinitionView._dictionaryDaya =  data;

        if (_.has(data, 'title')) {
          _tableDefinitionView._breadcrumbName = data.title;

          _tableDefinitionView._prettyName = (_DICTIONARY_CONSTANTS.DICTIONARY_ENTITY_MAP[data.category] || 'Unknown') +
                                             ': ' + data.title + ' Definition';


        }
        console.log(_tableDefinitionView._prettyName);

        return _tableDefinitionView;
      };

      _tableDefinitionView.setDictionaryData(_tableDefinitionView._dictionaryData);
    }

    // Inherit our defaults for View
    TableDefinitionsView.prototype = new View();

    TableDefinitionsView.prototype.renderDefinitionView = function(currentDictionary) {
      var _tableDefinitionView = this;
      console.log(currentDictionary);
    };

    return TableDefinitionsView;

  })();



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
        _ID: 'TABLE',
        ENTITY_LIST: 'table-entity-list',
        TERM_DEFINITION: 'table-definition-view'
      }
    },
    VIEW_STATE: {
      ENTER: 'enter', EXIT: 'exit', RENDERED: 'rendered'
    },
    VIEW_UPDATE_EVENT_TYPES: {
      INIT: 'initialize-view',
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
  var _urlParamsCache = null;

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

    if (! viewObject instanceof View) {
      console.warn('View object is not an instance of View - something bad might happen moving forward!');
    }

    _viewObject.view = viewObject;

    _viewObject.eventType = viewEventType || _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.DEFAULT;
    _viewObject.params = viewParams || null;
  }

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
      view: new TableEntityListView(tableViews.summary, dictionaryData, actionCallbackFn)
    };

    if (_.has(urlParams, 'id')) {
      dictionaryData = dictionaryData.dictionaryMap[urlParams.id] || null;
    }

    views[_DICTIONARY_CONSTANTS.VIEWS.TABLE._ID][_DICTIONARY_CONSTANTS.VIEWS.TABLE.TERM_DEFINITION] = {
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

    if (_urlParamsCache) {
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