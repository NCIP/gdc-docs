$(function () {
    var _searchItemClass = 'hp-search__item';
    var _VALID_QUERY_LENGTH = 3;
    var $cancelButton = $('.hp-search__cancel');
    var $inputBox = $('.hp-search__input');
    var $results = $('.hp-search__results');
    var $resultsWrapper = $('.hp-search__wrapper-results');
    var $resultsContainer = $('.hp-search__results-container');
    var $searchContainer = $('.hp-search');
    var $searchContentBody = $('.hp-search__body');
  
    $inputBox.focus();
  
    function _debounce(func, wait, immediate) {
      var _timeout;
      return function () {
        var context = this;
        var args = arguments;
        var later = function () {
          _timeout = null;
          if (!immediate) func.apply(context, args);
        };
        var callNow = immediate && !_timeout;
        clearTimeout(_timeout);
        _timeout = setTimeout(later, wait);
        if (callNow) {
          func.apply(context, args);
        }
      };
    }
  
    $.get(base_url + '/mkdocs/search_index.json', function (data) {
      var index = lunr(function () {
        this.field('title', {
          boost: 10
        });
        this.field('text');
        this.ref('location');
      });
  
      var documents = {};
      var doc;
  
      for (var i = 0; i < data.docs.length; i++) {
        doc = data.docs[i];
        doc.location = base_url + doc.location;
        index.add(doc);
        documents[doc.location] = doc;
      }
  
      function _search() {
        var query = $.trim($inputBox.val());
        $results.empty();
  
        if (query.length > 0) {
          $cancelButton.show();
        } else {
          $cancelButton.hide();
        }
  
        if (query.length < _VALID_QUERY_LENGTH || query === '') {
          $resultsWrapper.hide();
          $searchContainer.removeClass('search-active');
          return;
        }
  
        $resultsWrapper.show();
        $searchContainer.addClass('search-active');
  
        var results = index.search(query);
        var resultsHTML = '';
  
        $searchContentBody.html('<strong><i class="fa fa-file-o"></i> ' + results.length + '</strong> ' + (results.length === 1 ? 'result' : 'results') + ' found for <strong>' + query + '</strong>');
  
        if (results.length === 0) {
          $resultsContainer.hide();
        } else {
          $resultsContainer.show();
          for (var i = 0; i < results.length; i++) {
            var result = results[i];
            var resultDoc = documents[result.ref];
            resultDoc.base_url = base_url;
            resultDoc.summary = resultDoc.text.substring(0, 200);
  
            resultsHTML += '<a class="' + _searchItemClass + ' animated fadeInTop" href="' + resultDoc.location + '">' +
              '<div class="search-body">' +
              '<h3>' + resultDoc.title + '</h3>' +
              '<p>' + resultDoc.summary + '</p>' +
              '</div></a>';
          }
  
          $results.append(resultsHTML);
          $results.highlight(query);
  
          setTimeout(function () { $('.' + _searchItemClass).removeClass('animated fadeInTop'); }, 500);
        }
      }
  
      $inputBox.on('keyup', _debounce(_search, 300));
  
      $cancelButton.on('click', function () {
        $inputBox.val('');
        _search();
      })
    });
  
  });