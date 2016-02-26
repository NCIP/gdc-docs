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

      _tableDefinitionView.setDictionaryData = function(data) {

        _tableDefinitionView._dictionaryData =  data;

        if (_.has(data, 'title')) {
          _tableDefinitionView._breadcrumbName = data.title;
          _tableDefinitionView._prettyName =  data.title;
        }

        return _tableDefinitionView;
      };

      _tableDefinitionView.setDictionaryData(_tableDefinitionView._dictionaryData);
    }

    ////////////////////////////////////////////////////////////
    // Custom Render
    ////////////////////////////////////////////////////////////
    TableDefinitionsView.prototype.renderView = function() {
      var _tableDefinitionView = this;

      if (_tableDefinitionView._dictionaryData) {
        _tableDefinitionView.renderDefinitionView(_tableDefinitionView._dictionaryData);
      }

      console.log('TableDefinitionsView Rendering!');

    };

    TableDefinitionsView.prototype.renderDefinitionView = function(currentDictionary) {
      var _tableDefinitionView = this;

      console.log(currentDictionary);

      _tableDefinitionView.renderHeader();
      _tableDefinitionView.renderSummaryTable();
      _tableDefinitionView.renderLinksTable();
      _tableDefinitionView.renderPropertiesTable();
    };


    TableDefinitionsView.prototype.renderHeader = function() {
      var _tableDefinitionView = this;
      _tableDefinitionView._d3ContainerSelection.append('h1').text(_tableDefinitionView.getPrettyName());
    };

    TableDefinitionsView.prototype.renderSummaryTable = function() {
      var _tableDefinitionView = this;

      var summaryTableContainerSel = _tableDefinitionView._d3ContainerSelection.append('div')
          .classed('dictionary-summary-table-container dictionary-definition-container', true);

      _renderSummaryTable(_tableDefinitionView, summaryTableContainerSel);

    };

    TableDefinitionsView.prototype.renderLinksTable = function() {
      var _tableDefinitionView = this;

      var linksTableContainerSel = _tableDefinitionView._d3ContainerSelection.append('div')
        .classed('dictionary-links-table-container dictionary-definition-container', true);

      _renderLinksTable(_tableDefinitionView, linksTableContainerSel);
    };

    TableDefinitionsView.prototype.renderPropertiesTable = function() {
      var _tableDefinitionView = this;

      var propertiesTableContainerSel = _tableDefinitionView._d3ContainerSelection.append('div')
        .classed('dictionary-properties-table-container dictionary-definition-container', true);

      _renderPropertiesTable(_tableDefinitionView, propertiesTableContainerSel);
    };

    ///////////////////////////////////////////////////////////
    // Render Helpers
    ///////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////
    // Properties Table
    ///////////////////////////////////////////////////////////////////////////////////////
    function _getPropertyValueRecursive(propertyVal) {

      //console.log(propertyVal);

      if (_.has(propertyVal, 'enum')) {
        return propertyVal.enum;
      }
      else if (_.has(propertyVal, 'items.properties')) {
        return _.keys(propertyVal.items.properties)
      }
      else if (_.has(propertyVal, 'properties')) {
        return _.keys(propertyVal.properties)
      }
      else if (_.has(propertyVal, 'type')) {
        if (_.has(propertyVal, 'format')) {
          return  propertyVal.type + ' (Format: ' + propertyVal.format + ')';
        }
        return propertyVal.type;
      }

      return propertyVal;
    }

    function _getPropertyValueOrType(property) {
      var typesValuesPrecedence = ['enum', 'anyOf', 'oneOf', 'type'],
          value = null;

      if (! _.isObject(property)) {
        return property;
      }

      for (var i = 0; i < typesValuesPrecedence.length; i++) {
        var propertyName = typesValuesPrecedence[i],
            propertyVal = property[propertyName];

        if (_.has(property, propertyName)) {
          var normalizedPropertyName = propertyName.toLowerCase();
          switch (normalizedPropertyName) {
            case 'enum':
            case 'type':
              value = {propertyName: normalizedPropertyName === 'enum' ? 'Enumeration' : propertyName, propertyValue: propertyVal};
              break;
            case 'anyof':
            case 'oneof':
              value = {propertyName: normalizedPropertyName === 'oneof' ? 'One of' : 'Any of', propertyValue: []};

              for (var j = 0; j < propertyVal.length; j++) {
                value.propertyValue.push(_getPropertyValueRecursive(propertyVal[j]));
              }
              break;
          }

          break;
        }
      }

      return value;
    }

    function _prepPropertiesTableData(dictionaryData) {
      var propertyData = [];

      var dictionaryProperties = _.get(dictionaryData, 'properties', false),
          requiredProperties = _.get(dictionaryData, 'required', false);




      var propertyIDs;

      var requiredPartition = _.partition(_.keys(dictionaryProperties), function(propertyName) {

        return requiredProperties && requiredProperties.indexOf(propertyName) >= 0;
      });

      if (requiredPartition.length === 2) {
        propertyIDs = _.sortBy(requiredPartition[0]).concat(_.sortBy(requiredPartition[1]));
      }
      else {
        propertyIDs = requiredPartition;
      }

      console.log(requiredPartition);

      for (var i = 0; i < propertyIDs.length; i++) {
        var p = [],
            propertyName = propertyIDs[i],
            property = dictionaryProperties[propertyName],
            description = property.description,
            valueOrType = _getPropertyValueOrType(property),
            CDE = _.get(dictionaryData, 'terms.' + propertyName + '.termDef'),
            isRequired = requiredProperties && requiredProperties.indexOf(propertyName) >= 0 ? 'Yes' : 'No';

        // Ignore system properties for now...
        if (dictionaryData.systemProperties.indexOf(propertyName) >= 0) {
         console.log('Skipping system property: ' + propertyName);
          continue;
        }


        p.push(propertyName);
        p.push(_valueOrDefault(description));
        p.push(_valueOrDefault(valueOrType));
        p.push(isRequired);
        p.push(_valueOrDefault(CDE)) ;
        propertyData.push(p);
      }

      //console.log(propertyValues);

      console.log('Property: ', propertyData);

      if (propertyData.length === 0) {
        return [_.times(5, _.constant(_DICTIONARY_CONSTANTS.DATA_FORMATS.MISSING_VAL))];
      }

      return propertyData;
    }

    function _renderPropertiesTable(_tableDefinitionView, tableContainerSelection) {
      var dictionaryData = _tableDefinitionView._dictionaryData;

      console.log(dictionaryData);

      tableContainerSelection.append('h2')
        .append('a')
        .attr('id', 'properties-table')
        .attr('href',  '#?view=' + _tableDefinitionView._name + '&id=' + dictionaryData.id + '&anchor=properties-table')
        .text('Properties');

      var definitionTable = tableContainerSelection.append('table')
        .classed('dictionary-properties-table', true);

      var tHead = definitionTable.append('thead'),
          tBody = definitionTable.append('tbody');

      tHead.append('tr')
        .classed('dictionary-properties-header', true)
        .selectAll('th')
        .data(['Property', 'Description', 'Acceptable Types or Values', 'Required?', 'CDE'])
        .enter()
        .append('th')
        .text(function(d) { return d; });

      var dataRows = _prepPropertiesTableData(dictionaryData);

      var tRows = tBody.selectAll('tr')
        .data(dataRows)
        .enter()
        .append('tr');

      tRows.selectAll('td')
        .data(function(row) {
          return row;
        })
        .enter()
        .append('td')
        .classed('required-val',function(d, i) {
          if (i === 3) {
            return d === 'Yes';
          }

          return false;
        })
        .html(function(d, i) {

          var data = d;

          if (i === 0 && _.isString(data)) {
            data = '<div id="' + data + '"><a  href="#?view=' + _tableDefinitionView.getViewName() + '&id='+ dictionaryData.id + '&anchor=' + data + '"><i class="fa fa-gear"></i> ' + data + '</a></div>';
          }

          if (_.isString(data)) {
            return data;
          }

          if (i === 4) {
            var cdeStr = '';

            if (data.term_url) {
              cdeStr = '<a href="' + (data.term_url ? data.term_url : '#') + '" target="_blank">' + data.cde_id +
                       '</a>' + ' - ' + data.source;
            }
            else {
              cdeStr = _valueOrDefault();
            }

            return cdeStr;
          }

          console.log(data);

          if (data.propertyName === 'type') {

            if (!_.isArray(data.propertyValue)) {
              return data.propertyValue;
            }
            else {
              return data.propertyValue.join(', ');
            }
          }

          var bullets = _.map(data.propertyValue, function (val, i) {

            var arrayVal = '';

            //console.log(val, i);

            if (_.isArray(val) && val.length === 1) {
              arrayVal = val[0];
            }
            else if (_.isArray(val) && val.length > 1) {
              arrayVal = val.join(', ');
            }
            else {
              arrayVal = val;
            }

            return '<li>' + arrayVal + '</li>';
          }).join('\n\t');


          data = '<ul class="bullets"><li>' + data.propertyName + ': <ul>' + bullets + '</ul></li></ul>';

          return data;
        });

  }

    function _renderSummaryTable(_tableDefinitionView, tableContainerSelection) {
      var dictionaryData = _tableDefinitionView._dictionaryData,
        category = _.get(_DICTIONARY_CONSTANTS.DICTIONARY_ENTITY_MAP, dictionaryData.category.toLowerCase(), dictionaryData.category),
        uniqueKeys = _.get(dictionaryData, 'uniqueKeys', [_DICTIONARY_CONSTANTS.DATA_FORMATS.MISSING_VAL]);

      tableContainerSelection.append('h2')
        .append('a')
        .attr('id', 'summary-table')
        .attr('href', '#?view=' + _tableDefinitionView._name + '&id=' + dictionaryData.id + '&anchor=summary-table')
        .text('Summary');

      var definitionTable = tableContainerSelection.append('table')
        .classed('dictionary-summary-table', true);

      var tHead = definitionTable.append('thead'),
        tBody = definitionTable.append('tbody');

      tHead.append('tr')
        .classed('dictionary-summary-header', true)
        .selectAll('th')
        .data(['Title', _tableDefinitionView.getPrettyName()])
        .enter()
        .append('th')
        .text(function (d) {
          return d;
        });


      var dataRows = [
        {id: 'category', title: 'Category', value: category},
        {id: 'description', title: 'Description', value: dictionaryData.description},
        {id: 'keys', title: 'Unique Keys', value: uniqueKeys}
      ];

      var tRows = tBody.selectAll('tr')
        .data(dataRows)
        .enter()
        .append('tr');

      tRows.selectAll('td')
        .data(function (row) {
          return [{id: row.id, value: row.title}, {id: row.id, value: row.value}];
        })
        .enter()
        .append('td')
        .html(function (d, i) {

          var data = d.value;

          if (i !== 1 || d.id !== 'keys') {
            return data;
          }

          if (_.isArray(data)) {

            var newData = '<ul class="bullets">' +
                          _.map(data, function (val) {
                            var arrayVal = '';

                            if (_.isArray(val) && val.length === 1) {
                              arrayVal = val[0];
                            }
                            else if (_.isArray(val) && val.length > 1) {
                              arrayVal = val.join(', ');
                            }
                            else {
                              arrayVal = val;
                            }

                            return '<li>' + arrayVal + '</li>';
                          }).join('\n') +
                          '</ul>';

            data = newData;

          }

          return data;
        });
    }


    ///////////////////////////////////////////////////////////////////////////////////////
    // Links Table
    ///////////////////////////////////////////////////////////////////////////////////////
    function createLinkData(link) {
      var linkData = [];

      if (! _.has(link, 'name')) {
        return _.times(3, _.constant(_DICTIONARY_CONSTANTS.DATA_FORMATS.MISSING_VAL));
      }

      var linkID = _.get(link, 'target_type'),
          linkLabel = link.label ?  link.label.split('_').join(' ') : _valueOrDefault();

      linkData.push({id: linkID, name: link.name});
      linkData.push(_capitalizeWords(_valueOrDefault(link.backref).split('_').join(' ')) + ' ' +
                    '<strong>' +   _capitalizeWords(linkLabel) + '</strong> '  +
                    _capitalizeWords(_valueOrDefault(link.target_type.split('_').join(' '))));
      linkData.push(link.required === true ? 'Yes' : 'No');

      return linkData;
    }


    function _prepLinksTableData(dictionaryData) {
      var transformedData = [];

      // Create Table Row Data
      var links = _.get(dictionaryData, 'links', false);

      if (! links || ! _.isArray(links) || links.length === 0) {
        return [createLinkData()];
      }

      for (var i = 0; i < links.length; i++) {
        var link = links[i];

        var linkSubgroups = _.get(link, 'subgroup', []);

        if (linkSubgroups.length > 0) {
          var subLinkData = [[],[],[]];

          for (var j = 0; j < linkSubgroups.length; j++) {
            var subLinkDataNode = createLinkData(linkSubgroups[j]);

            if (_.isArray(subLinkDataNode)) {
              subLinkData[0].push(subLinkDataNode[0]);
              subLinkData[1].push(subLinkDataNode[1]);
              // Link dictates whether subgroup is required
              subLinkData[2].push(link.required === true ? 'Yes': 'No');
            }
          }

          if (subLinkData.length > 0) {
            transformedData.push(subLinkData);
          }

          continue;
        }


        var linkDataNode = createLinkData(link);

        if (linkDataNode.length > 0) {
          transformedData.push(linkDataNode);
        }

      }

      return transformedData;
    }

    function _renderLinksTable(_tableDefinitionView, tableContainerSelection) {
      var dictionaryData = _tableDefinitionView._dictionaryData;

      tableContainerSelection.append('h2')
        .append('a')
        .attr('id', 'links-table')
        .attr('href',  '#?view=' + _tableDefinitionView._name  + '&id=' + dictionaryData.id + '&anchor=links-table')
        .text('Links');

      var definitionTable = tableContainerSelection.append('table')
        .classed('dictionary-links-table', true);

      var tHead = definitionTable.append('thead'),
        tBody = definitionTable.append('tbody');

      tHead.append('tr')
        .classed('dictionary-links-header', true)
        .selectAll('th')
        .data(['Links to', 'Relationship', 'Required?'])
        .enter()
        .append('th')
        .text(function(d) { return d; });


      var dataRows = _prepLinksTableData(dictionaryData);

      var tRows = tBody.selectAll('tr')
        .data(dataRows)
        .enter()
        .append('tr');

      tRows.selectAll('td')
        .data(function(row) {
          return row;
        })
        .enter()
        .append('td')
        .classed('required-val',function(d, i) {
          if (i === 2) {
            if (_.isArray(d)) {
              return _.first(d) === 'Yes'
            }

            return d === 'Yes';
          }

          return false;
        })
        .html(function(data, i) {

          if (i === 0) {

            var link = data;

            if (_.isString(link)) {
              return link;
            }

            var isNotSubgroup = false;

            if (!_.isArray(link)) {
              isNotSubgroup = true;
              link = [data];
            }

            link = _.map(link, function(l) {
              return  '<a href="#?view=' + _tableDefinitionView.getViewName() + '&id=' + l.id + '&_top=1" title="' +
                      (isNotSubgroup ? 'Entity' : 'Entity Subgroup') + '">' +
                      (isNotSubgroup ? '<i class="fa fa-file-o"></i>': '<i class="fa fa-sitemap"></i>') + ' ' +
                      _capitalizeWords(_capitalizeWords(l.id.split('_').join(' '))) +
                      '</a>';
            }).join('<br />\n');


            return link;
          }


          if (_.isString(data)) {
            return data;
          }

          if (i === 2) {
            return _.first(data);
          }


          var newData = data.join('<br />');

          //console.log(newData);

          return newData;
        });

    }



    return TableDefinitionsView;

  })();



  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////


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

    }

    TableEntityListView.prototype.renderView = function () {
      var _tableEntityListView = this;
      var categoryMap = _tableEntityListView._dictionaryData.dictionaryMapByCategory;
      var categoryKeys = _DICTIONARY_CONSTANTS.ENTITY_LIST_DICTIONARY_KEY_ORDER;
      console.log(categoryKeys);

      for (var i = 0; i < categoryKeys.length; i++) {
        var category = categoryKeys[i];

        _tableEntityListView.renderEntity(category, categoryMap[category]);

      }

      console.log('TableEntityListView Rendering!');


      _tableEntityListView._state = _DICTIONARY_CONSTANTS.VIEW_STATE.RENDERED;
      _tableEntityListView._callbackFn.call(null, new Dictionary._ViewUpdateObject(_tableEntityListView));
    };


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
          case 'case':
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

          // Here we just want to fire the event programmatically
          if ( _DICTIONARY_CONSTANTS.BROWSER_CAPABILITIES.HASH_CHANGE_EVENT ) {
            return;
          }

          _tableEntityListView._callbackFn.call(
            null, new Dictionary._ViewUpdateObject(_tableEntityListView,  _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.INNER_NAV, {id: category})
          );
        })
        .html(function() {
          var tooltipText = getTooltipText();
          return '<i class="fa fa-book"></i> ' + _.get(_DICTIONARY_CONSTANTS.DICTIONARY_ENTITY_MAP, category.toLowerCase(), category) +
                 (_.isString(tooltipText) ? '<span><i></i>' + tooltipText + '</span> &nbsp;<i style="color: #ccc;" class="fa fa-info-circle"></i>' : '');
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

          // Don't fire event twice since clicking on an href creates a pop event
          if ( _DICTIONARY_CONSTANTS.BROWSER_CAPABILITIES.HASH_CHANGE_EVENT ) {
            return;
          }

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



  function _valueOrDefault(val) {
    return val !== null && typeof val !== 'undefined' ? val : _DICTIONARY_CONSTANTS.DATA_FORMATS.MISSING_VAL;
  }

  function _capitalizeWords(str) {
    return str.replace(/[^\s]+/g, function(word) {
      return word.replace(/^[a-z]/i, function(firstLetter) {
        return firstLetter.toUpperCase();
      });
    });
  }

  /////////////////////////////////////////////////////////
  // Parent View Definition
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

    if (typeof _view.renderView === 'undefined') {
      throw Error('You must define your own renderView method in your view!');
    }

    /////////////////////////////////////////////////////////
    // Public View API
    /////////////////////////////////////////////////////////
    _view.render = function() {
      console.log('Rendering View!');

      _view._d3ContainerSelection.html('');

      // Template method - inherited functions define this!
      _view.renderView();

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
    };

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




})(window.Dictionary._DICTIONARY_CONSTANTS, window.Dictionary._);