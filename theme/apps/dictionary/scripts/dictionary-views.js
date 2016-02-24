(function(_DICTIONARY_CONSTANTS, _) {
  'use strict';

  if (typeof window.Dictionary !== 'function') {
    console.warn('Could not find the Dictionary Global - are you sure you included the main dictionary JS File.\nAborting View Init.');
    return;
  }

  Dictionary._Views = {};

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
      _view._callbackFn.call(null, new Dictionary._ViewUpdateObject(this));

      return _view;
    };

    _view.show = function() {
      console.log('Showing!');

      _view._isHidden = false;
      _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.ENTER;
      _view._d3ContainerSelection.style('display', 'block').transition().duration(250).style('opacity', 1);

      _view._callbackFn.call(null, new Dictionary._ViewUpdateObject(this));

      return _view;
    };

    _view.hide = function() {
      console.log('Hiding!');

      _view._isHidden = true;
      _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.EXIT;
      _view._d3ContainerSelection.transition().duration(250).style('opacity', 0);
      setTimeout(function() { _view._d3ContainerSelection.style('display', 'none'); }, 250);


      _view._callbackFn.call(null, new Dictionary._ViewUpdateObject(this));

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

  Dictionary._ViewUpdateObject = function(viewObject, viewEventType, viewParams) {
    var _viewObject = this;

    _viewObject.view = viewObject;

    _viewObject.eventType = viewEventType || _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.DEFAULT;
    _viewObject.params = viewParams || null;
  };


  /////////////////////////////////////////////////////////
  // TableEntityListView
  /////////////////////////////////////////////////////////
  Dictionary._Views.TableEntityListView = (function() {


    function TableEntityListView() {

      var _tableEntityListView = this;
      // Inherit from View
      View.apply(_tableEntityListView, arguments);

      _tableEntityListView._parentViewName = _DICTIONARY_CONSTANTS.VIEWS.TABLE._ID;
      _tableEntityListView._name = _DICTIONARY_CONSTANTS.VIEWS[_tableEntityListView._parentViewName].ENTITY_LIST;
      _tableEntityListView._prettyName = 'Entity List';

      /////////////////////////////////////////
      // Custom Render
      /////////////////////////////////////////
      _tableEntityListView.render = function () {

        var categoryMap = _tableEntityListView._dictionaryData.dictionaryMapByCategory;

        for (var category in categoryMap) {
          if (categoryMap.hasOwnProperty(category)) {
            _tableEntityListView.renderEntity(category, categoryMap[category]);
          }
        }

        console.log('TableEntityListView Rendering!');


        _tableEntityListView._state = _DICTIONARY_CONSTANTS.VIEW_STATE.RENDERED;
        _tableEntityListView._callbackFn.call(null, new Dictionary._ViewUpdateObject(_tableEntityListView));
      };

    }


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
            null, new Dictionary._ViewUpdateObject(_tableEntityListView,  _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INNER_NAV, {id: category})
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
            null, new Dictionary._ViewUpdateObject(_tableEntityListView,  _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.NAV, {id: data.id})
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
  Dictionary._Views.TableDefinitionsView = (function() {
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
        _tableDefinitionView._callbackFn.call(null, new Dictionary._ViewUpdateObject(_tableDefinitionView));
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


    TableDefinitionsView.prototype.renderDefinitionView = function(currentDictionary) {
      var _tableDefinitionView = this;
      console.log(currentDictionary);
    };

    return TableDefinitionsView;

  })();



  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////




})(window.Dictionary._DICTIONARY_CONSTANTS, window.Dictionary._);