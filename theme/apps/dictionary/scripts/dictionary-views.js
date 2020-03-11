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
      var headerSelection = _tableDefinitionView._d3ContainerSelection.append('h1')
            .text(_tableDefinitionView.getPrettyName());

      var definitionControlsSelection = headerSelection.append('div')
            .classed('definition-controls-container', true);


      // Exclude the below from download
      var excludeCategories = _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_DOWNLOAD_BLACKLIST;

      if (excludeCategories.indexOf(_tableDefinitionView._dictionaryData.ui_category.toLowerCase()) < 0) {

        var updateHREFFunction = function() {
          d3.select(this)
            .attr('href', _tableDefinitionView._parentDictionary.getDictionaryTemplateURL(_tableDefinitionView._dictionaryData.id));
        };

        definitionControlsSelection
          .append('a')
          .attr('href', '')
          .on('mouseenter', updateHREFFunction)
          .on('focus', updateHREFFunction)
          .attr('title', 'Download the ' + _tableDefinitionView.getPrettyName() + ' template.')
          .classed('dictionary-control-bttn dictionary-template-download-bttn', true)
          .html('<span aria-hidden="true" class="fa fa-cloud-download"></span> &nbsp;Download Template');


        _tableDefinitionView._parentDictionary
          .fetchDictionaryTemplate('views/entity-definition-controls.view.html')
          .then(function(html) {
            var templateDefinitionControlSelection = definitionControlsSelection.append('div')
              .style({display: 'inline-block', 'margin-left': '2rem'})
              .html(html);

            var viewDataFormatLabel = templateDefinitionControlSelection.select('.data-format-value');

            // Set the label to the current default on control load
            viewDataFormatLabel.text(_tableDefinitionView._parentDictionary.getDefaultDictionaryTemplateDownloadFormat().toUpperCase());


            templateDefinitionControlSelection.selectAll('a').on('click', function() {
              var option = d3.select(this);

              var dataFormat = option.attr('data-format') || _DICTIONARY_CONSTANTS.ENDPOINT_PARAMS.TSV_FORMAT,
                  dataFormatLabelVal = dataFormat.toUpperCase();

              viewDataFormatLabel.text(dataFormatLabelVal);

              _tableDefinitionView._parentDictionary.setDefaultDictionaryTemplateDownloadFormat(dataFormat);
            });

          });

    }

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
              value = {
                propertyName: normalizedPropertyName === 'enum' ? 'Enumeration' :
                  (propertyVal === 'boolean' ? 'boolean' : propertyName ),
                propertyValue: propertyVal === 'boolean' ? ['true', 'false'] : propertyVal
              };
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
          requiredProperties = _.get(dictionaryData, 'required', false),
          // Exclude System Properties
          excludeProperties = _.get(dictionaryData, 'systemProperties', []);

      // Exclude the type property
      excludeProperties = excludeProperties.concat(_DICTIONARY_CONSTANTS.PROPERTY_EXCLUDES);

      // Exclude unique keys
      if (_.isArray(dictionaryData.uniqueKeys)) {

        var uniqueKeys = _.reduce(dictionaryData.uniqueKeys, function(flattenedKeys, keys) {
          flattenedKeys = flattenedKeys.concat(keys);
          return flattenedKeys;
        }, []);

        excludeProperties = excludeProperties.concat(uniqueKeys);
      }


      var links = _prepLinksTableData(dictionaryData);

      // Exclude links in properties
      if (links !== null) {

        if (links.topLevelLinks.length) {
          excludeProperties = excludeProperties.concat(_.map(links.topLevelLinks, function(l) { return l[0].name; }));
        }

        // Exclude Sublinks in properties
        if (links.subLinks.length) {
          var sublinks = links.subLinks;

          var justSublinkNames = _.map(sublinks, function(sublinkIDObjs) {
            // Link object in the first index of each sublink array group
            return _.reduce(_.first(sublinkIDObjs), function(sublinksArray, l) {
              sublinksArray.push(l.name);
              return sublinksArray;
            }, []);

          });

          var flattenedSublinkIDs = _.reduce(justSublinkNames, function(propertyList, names) {
            propertyList = propertyList.concat(names);
            return propertyList;
          }, []);

          excludeProperties = excludeProperties.concat(flattenedSublinkIDs);
        }

        if (links.excludedLinks.length) {
          var excludedLinks = links.excludedLinks;
          var flattenedExcludedLinks = _.reduce(excludedLinks, function(list, link){
            list.push(link.name);
            return list;
          }, []);

          excludeProperties = excludeProperties.concat(flattenedExcludedLinks);
        }
      }


      // Make the excluded properties unique
      excludeProperties = _.uniq(excludeProperties);

      console.log('Excluded Properties: ', excludeProperties);

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

      var dictionaryTerms = dictionaryData.terms || {};

      for (var i = 0; i < propertyIDs.length; i++) {
        var p = [],
            propertyName = propertyIDs[i],
            property = dictionaryProperties[propertyName],
            description = _.get(property, 'term.description', false) || _.get(property, 'description', false) || _.get(dictionaryTerms, propertyName + '.description', 'n/a'),
            valueOrType = _getPropertyValueOrType(property),
            CDE = _.get(dictionaryTerms, propertyName + '.node.property.termDef') ||
              _.get(dictionaryTerms, propertyName + '.common.termDef'),
            isRequired = requiredProperties && requiredProperties.indexOf(propertyName) >= 0 ? 'Yes' : 'No';

        // Ignore system properties for now...
        if (excludeProperties.indexOf(propertyName) >= 0) {
          console.log('Skipping excluded property: ' + propertyName);
          continue;
        }

        p.push(propertyName);
        p.push(_valueOrDefault(description));
        p.push(_valueOrDefault(valueOrType));
        p.push(isRequired);
        p.push(_valueOrDefault(CDE)) ;
        propertyData.push(p);
      }

      if (propertyData.length === 0) {
        return null; //[_.times(5, _.constant(_DICTIONARY_CONSTANTS.DATA_FORMATS.MISSING_VAL))];
      }

      return propertyData;
    }

    function _renderPropertiesTable(_tableDefinitionView, tableContainerSelection) {
      var dictionaryData = _tableDefinitionView._dictionaryData;

      tableContainerSelection.append('h2')
        .append('a')
        .attr('id', 'properties-table')
        .attr('href',  '#?view=' + _tableDefinitionView._name + '&id=' + dictionaryData.id + '&anchor=properties-table')
        .text('Properties');

      var dataRows = _prepPropertiesTableData(dictionaryData);

      if (dataRows === null) {
        tableContainerSelection.append('h3')
          .style({'margin-left': '1rem', 'color': '#777'})
          .html('<i class="fa fa-gears"></i> &nbsp;No Properties for this entity.');

        return;
      }

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
        .classed('no-print',function(d, i) {
          if (i === 4) {
            return true;
          }

          return false;
        })
        .text(function(d) { return d; });

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
        .classed('no-print',function(d, i) {
          if (i === 4) {
            return true;
          }

          return false;
        })
        .html(function(d, i) {

          var data = d;

          switch(i) {
            case 0:
              if (_.isString(data)) {

                if (data === _DICTIONARY_CONSTANTS.DATA_FORMATS.MISSING_VAL) {
                  return data;
                }

                data = '<div id="' + data + '"><a class="monospace dictionary-anchor"  href="#?view=' + _tableDefinitionView.getViewName() + '&id='+ dictionaryData.id + '&anchor=' + data + '"><i class="fa fa-gear"></i> ' + data + '</a></div>';
              }
              break;
            case 1:
              return '<div class="property-description">' + data + '</div>';
              break;
            case 4:
              var cdeStr = '';

              if (data.term_url) {
                cdeStr = '<a href="' + (data.term_url ? data.term_url : '#') + '" target="_blank">' + data.cde_id +
                         '</a>' + ' - ' + data.source;
              }
              else {
                cdeStr = _valueOrDefault();
              }

              return cdeStr;
              break;
            default:
              break;

          }


          if (_.isString(data)) {
            return data;
          }

          if (data.propertyName === 'type') {

            if (!_.isArray(data.propertyValue)) {
              return data.propertyValue;
            }
            else {
              return data.propertyValue.join(', ');
            }
          }

          var bullets = _.map(data.propertyValue, function (val) {

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
          }).join('\n\t');


          data = '<div class="values-accordion"><ul class="bullets monospace"><li>' + data.propertyName + ': <ul>' + bullets + '</ul></li></ul></div>';

          return data;
        });

  }

    function _renderSummaryTable(_tableDefinitionView, tableContainerSelection) {
      var dictionaryData = _tableDefinitionView._dictionaryData,
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

      var dataRows = [
        {id: 'type', title: 'Type', value: dictionaryData.id},
        {id: 'category', title: 'Category', value: dictionaryData.category},
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

          if (i === 1 && d.id === 'type') {
            return '<span class="monospace">' + data + '</span>';
          }

          if (i !== 1 || d.id !== 'keys') {
            return data;
          }

          if (_.isArray(data)) {

            var newData = '<ul class="bullets monospace">' +
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
        return _.times(4, _.constant(_DICTIONARY_CONSTANTS.DATA_FORMATS.MISSING_VAL));
      }

      var linkID = _.get(link, 'target_type'),
          linkLabel = link.label ?  link.label.split('_').join(' ') : _valueOrDefault();

      linkData.push({id: linkID, name: link.name});
      linkData.push(link.name);
      linkData.push(_capitalizeWords(_valueOrDefault(link.backref).split('_').join(' ')) + ' ' +
                    '<strong>' +   _capitalizeWords(linkLabel) + '</strong> '  +
                    _capitalizeWords(_valueOrDefault(link.target_type.split('_').join(' '))));
      linkData.push(link.required === true ? 'Yes' : 'No');

      return linkData;
    }


    function _prepLinksTableData(dictionaryData) {
      var transformedData = [],
          topLevelLinks = [],
          subLinks = [],
          excludedLinks = [];

      // Create Table Row Data
      var links = _.get(dictionaryData, 'links', false);

      if (! links || ! _.isArray(links) || links.length === 0) {
        return null;
      }


      var exclusions = _DICTIONARY_CONSTANTS.LINK_EXCLUDES;

      for (var i = 0; i < links.length; i++) {
        var link = links[i],
            linkSubgroups = _.get(link, 'subgroup', []),
            isLinkIncluded = true;

        if (linkSubgroups.length > 0) {
          var subLinkData = [[],[],[],[]];

          // Sort sublinks by name
          var sortedLinkSubgroups = _.sortBy(linkSubgroups, function(l) {
            return l.name;
          });


          for (var j = 0; j < sortedLinkSubgroups.length; j++) {
            var subLinkDataNode = createLinkData(sortedLinkSubgroups[j]);

            isLinkIncluded = exclusions.indexOf(subLinkDataNode[0].id) < 0;

            if (_.isArray(subLinkDataNode) && isLinkIncluded) {
              subLinkData[0].push(subLinkDataNode[0]);
              subLinkData[1].push(subLinkDataNode[1]);
              subLinkData[2].push(subLinkDataNode[2]);
              // Link dictates whether subgroup is required
              subLinkData[3].push(link.required === true ? 'Yes': 'No');
            }

            if (! isLinkIncluded) {
              excludedLinks.push(subLinkDataNode[0]);
            }
          }

          if (subLinkData.length > 0) {
            subLinks.push(subLinkData);
          }

          continue;
        }


        var linkDataNode = createLinkData(link);

        isLinkIncluded = exclusions.indexOf(linkDataNode[0].id) < 0;

        if (! isLinkIncluded) {
          excludedLinks.push(linkDataNode[0]);
          continue;
        }

        if (linkDataNode.length > 0) {
          topLevelLinks.push(linkDataNode);
        }

      }

      if (topLevelLinks.length === 0 && subLinks.length === 0) {
        return {links: transformedData, topLevelLinks: topLevelLinks, subLinks: subLinks, excludedLinks: excludedLinks};
      }

      var sortASCandRequiredFirst = function (l) {
        return (l[2] === 'Yes' ? 'a' : 'z') + l[0].name;
      };

      if (topLevelLinks.length) {
        transformedData = _.sortBy(topLevelLinks, sortASCandRequiredFirst);
      }

      // If we have links and the first in group i
      if (transformedData.length && _.first(transformedData)[2] === 'Yes') {
        transformedData = transformedData.concat(subLinks);
      }
      else {
        transformedData = subLinks.concat(transformedData);
      }

      return {links: transformedData, topLevelLinks: topLevelLinks, subLinks: subLinks, excludedLinks: excludedLinks};
    }

    function _renderLinksTable(_tableDefinitionView, tableContainerSelection) {
      var dictionaryData = _tableDefinitionView._dictionaryData;

      tableContainerSelection.append('h2')
        .append('a')
        .attr('id', 'links-table')
        .attr('href',  '#?view=' + _tableDefinitionView._name  + '&id=' + dictionaryData.id + '&anchor=links-table')
        .text('Links');


      var dataRows = _prepLinksTableData(dictionaryData);

      if (dataRows === null || (dataRows.links.length === 0 && dataRows.subLinks.length === 0)) {
        tableContainerSelection.append('h3')
          .style({'margin-left': '1rem', 'color': '#777'})
          .html('<i class="fa fa-unlink"></i> &nbsp;No Links for this entity.');
        return;
      }

      dataRows = dataRows.links;

      var definitionTable = tableContainerSelection.append('table')
        .classed('dictionary-links-table', true);

      var tHead = definitionTable.append('thead'),
        tBody = definitionTable.append('tbody');

      tHead.append('tr')
        .classed('dictionary-links-header', true)
        .selectAll('th')
        .data(['Links to Entity', 'Link Name' , 'Relationship', 'Required?'])
        .enter()
        .append('th')
        .html(function(d, i) {
          if (i === 1) {
            return '<span class="dictionary-tooltip"><em>' + d  +  '</em><span class="dictionary-tooltip-content"><i></i>The links ' +
                   'should be included in the files uploaded to the GDC. ' +
                   'For more information, please refer to the ' +
                   '<a href="/Data_Submission_Portal/Users_Guide/Upload_Data/#step1-prepare-files" target="_blank">' +
                   '<!-- b class="fa fa-external-link"></b --> File ' +
                   'Preparation\'s User Guide</a>.</span></span>';
          }
          return d;

        });

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
            if (_.isArray(d)) {
              return _.first(d) === 'Yes'
            }

            return d === 'Yes';
          }

          return false;
        })
        .html(function(data, i) {

          switch (i) {
            case 0:
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
              break;
            case 1:

              if (_.isString(data)) {
                return data;
              }

              return '<span class="monospace">&middot; ' + data.join('<br />\n&middot; ') + '</span>';
              break;
            case 3:
              if (! _.isArray(data)) {
                return data;
              }

              return _.first(data);

              break;

            default:
              if (_.isString(data)) {
                return data;
              }

              var newData = data.join('<br />');

              return newData;
            break;
          }

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

      _tableEntityListView._parentDictionary
        .fetchDictionaryTemplate('views/entity-list-controls.view.html')
        .then(function(html) {
          var templateDefinitionControlSelection = d3.select(_DICTIONARY_CONSTANTS.VIEWS._STATIC.DICTIONARY_CONTROLS).append('div')
            .style({display: 'inline-block', 'margin-left': '2rem'})
            .html(html);

          var viewDataFormatLabel = templateDefinitionControlSelection.select('.data-format-value');

          // Set the label to the current default on control load
          viewDataFormatLabel.text(_tableEntityListView._parentDictionary.getDefaultDictionaryTemplateDownloadFormat().toUpperCase());


          templateDefinitionControlSelection.selectAll('a').on('click', function() {
            var option = d3.select(this);

            var dataFormat = option.attr('data-format') || _DICTIONARY_CONSTANTS.ENDPOINT_PARAMS.TSV_FORMAT,
                dataFormatLabelVal = dataFormat.toUpperCase();

            viewDataFormatLabel.text(dataFormatLabelVal);

            _tableEntityListView._parentDictionary.setDefaultDictionaryTemplateDownloadFormat(dataFormat);
          });

        });

      for (var i = 0; i < categoryKeys.length; i++) {
        var category = categoryKeys[i];

        _tableEntityListView.renderEntity(category, categoryMap[category], _getSortFnForCategory(category, categoryMap[category]));
      }

      // Print any remaining not explicitly sorted keys that may not be in the hardcoded order...
      var leftOverCategories = _.difference(
        _.keys(_tableEntityListView._dictionaryData.dictionaryMapByCategory), categoryKeys
      ).filter(function (key) {
        return _DICTIONARY_CONSTANTS.CATEGORY_EXCLUDES.indexOf(key) === -1
      });

      if (leftOverCategories.length) {
        console.warn('Sorted Category Differences: ', leftOverCategories);

        for (i = 0; i < leftOverCategories.length; i++) {
          var category = leftOverCategories[i];
          _tableEntityListView.renderEntity(category, categoryMap[category], _getSortFnForCategory(category, categoryMap[category]));
        }
      }

      console.log('TableEntityListView Rendering!');

      _tableEntityListView._state = _DICTIONARY_CONSTANTS.VIEW_STATE.RENDERED;
      _tableEntityListView._callbackFn.call(null, new Dictionary._ViewUpdateObject(_tableEntityListView));
    };


    TableEntityListView.prototype.renderEntity = function(category, categoryData, categorySortFn) {
      var _tableEntityListView = this;

      categoryData = categorySortFn.call(_tableEntityListView, categoryData);

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
          case 'administrative':
          case 'TBD':
            tooltipText = 'The entities listed in this category are maintained by GDC. Data submission is not applicable.';
            break;
          default:
            break;
        }

        return tooltipText;
      };

      var tHeadRow = tHead.append('tr')
        .append('th')
        .attr('colspan', 2)
        .classed('dictionary-entity-header', true);

      tHeadRow.append('a')
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
          return '<span class="fa fa-book" aria-hidden="true"></span> <em>' + _capitalizeWords(_.get(_DICTIONARY_CONSTANTS.DICTIONARY_ENTITY_MAP, category.toLowerCase(), category)) + '</em>' +
                 (_.isString(tooltipText) ? '<span class="dictionary-tooltip-content"><i></i>' + tooltipText + '</span> &nbsp;<!-- i style="color: #ccc;" class="fa fa-info-circle"></i -->' : '');
        });

      // Exclude the below from download
      var excludeCategories = _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_DOWNLOAD_BLACKLIST;

      if (excludeCategories.indexOf(category.toLowerCase()) < 0) {
        tHeadRow.append('div')
          .classed('dictionary-download-category-btn-container', true)
          .append('a')
          .attr('href', 'javascript:void(0)')
          .attr('title', 'Download All Templates for the ' +
                         _.get(_DICTIONARY_CONSTANTS.DICTIONARY_ENTITY_MAP, category.toLowerCase(), category) +
                         ' Category')
          .on('click', function () {


            // TODO: Remove hardcoding but for now this is necessary...
            var oneOffDictionaries = ['annotation', 'case'];

            if (oneOffDictionaries.indexOf(category) >= 0) {

              _tableEntityListView._callbackFn.call(
                null, new Dictionary._ViewUpdateObject(_tableEntityListView, _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.TEMPLATE_DOWNLOAD_REQUESTED, {
                  id: category
                })
              );

              return;
            }

            var exclusions = category === 'submittable_data_file' ?
              _tableEntityListView._dictionaryData.dictionaryMapByCategory.generated_data_file.map(function(f) { return f.id; }) :
              _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_EXCLUDES[category];
            var inclusions = _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_INCLUDES[category];

            _tableEntityListView._callbackFn(new Dictionary._ViewUpdateObject(
              _tableEntityListView,
              _DICTIONARY_CONSTANTS.VIEW_UPDATE_EVENT_TYPES.TEMPLATE_DOWNLOAD_BY_CATEGORY_REQUESTED,
              {
                id: category === 'submittable_data_file' ? 'data_file,metadata_file' : category,
                excludes: inclusions
                  ? _.difference(categoryData.map(function(x) { return x.id; }), inclusions) : exclusions
              })
            );
          })
          .html('<span aria-hidden="true" class="fa fa-cloud-download"></span>&nbsp; Download  &nbsp;');
      }

      var tRows = tBody.selectAll('tr')
        .data(
          categoryData.filter(function (entity) {
            return _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_INCLUDES[entity.category]
              ? _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_INCLUDES[entity.category].indexOf(entity.id) > -1
              : true;
          }).filter(function(d) { return _DICTIONARY_CONSTANTS.LINK_EXCLUDES.indexOf(d.id) === -1; })
          .sort(function (a, b) {
            return _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_INCLUDES[a.category]
              ? _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_INCLUDES[a.category].indexOf(a.id) -
                _DICTIONARY_CONSTANTS.CATEGORY_TEMPLATE_INCLUDES[b.category].indexOf(b.id)
              : 0;
          })
        )
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

  function _getSortFnForCategory(category, categoryData) {
    var sortFunction = _.constant(categoryData);

    switch(category.toLowerCase()) {
      case 'biospecimen':
        sortFunction = function(categoryData) {
          return _.orderBy(categoryData, ['title'], ['desc']);
        };
        break;
      default:
        break;
    }

    return sortFunction;
  }


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
  function View(d3ContainerSelection, dictionaryData, actionCallbackFn, parentDictionary) {
    var _view = this;

    _view._d3ContainerSelection = d3ContainerSelection;
    _view._dictionaryData = dictionaryData;
    _view._name = 'Unknown View';
    _view._callbackFn = actionCallbackFn || _.noop;
    _view._isHidden = true;
    _view._parentViewName = '';
    _view._prettyName = _view._name;
    _view._breadcrumbName = null;

    // Note: The below is optional and should only be used by views
    // that need access to the dictionary object's functionality
    _view._parentDictionary = parentDictionary || null;

    if (typeof _view.renderView === 'undefined') {
      throw Error('You must define your own renderView method in your view!');
    }

    /////////////////////////////////////////////////////////
    // Public View API
    /////////////////////////////////////////////////////////
    _view.render = function() {
      console.log('Rendering View!');

      _view._d3ContainerSelection.html('');

      // Show the controls before rendering
      d3.select(_DICTIONARY_CONSTANTS.VIEWS._STATIC.DICTIONARY_CONTROLS).html('');

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
      _view._d3ContainerSelection.style('display', 'block').classed('active', true).transition().duration(250).style('opacity', 1);

      _view._callbackFn.call(null, new Dictionary._ViewUpdateObject(this));

      return _view;
    };

    _view.hide = function() {
      console.log('Hiding!');

      _view._isHidden = true;
      _view._state = _DICTIONARY_CONSTANTS.VIEW_STATE.EXIT;
      _view._d3ContainerSelection.classed('active', false).transition().duration(250).style('opacity', 0);
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
