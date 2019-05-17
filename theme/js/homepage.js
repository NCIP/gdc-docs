$(function () {
  var $searchContainer = $('.hp-search');
  var $inputBox = $('.hp-search__input');
  var searchInput = document.getElementById('hp-search__input');
  var $resultsWrapper = $('.hp-search__wrapper-results');
  var $resultsContainer = $('.hp-search__results-container');
  var $results = $('.hp-search__results');
  var $searchContentBody = $('.hp-search__body');
  var $cancelButton = $('.hp-search__cancel');
  var _VALID_QUERY_LENGTH = 3;
  var _isSearchActive = false;
  var _searchItemClass = 'hp-search__item';

  $inputBox.focus();

  $inputBox.on('blur', function () {
    var $this = $(this);
    if ($this.val() === '') {
      $this.closest('.hp-search').removeClass('search-active');
    }
  });

  function _debounce(func, wait, immediate) {
    var _timeout;

    return function () {
      var context = this,
        args = arguments,
        later = function () {
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

  function _resetSearch() {
    $inputBox.val('');
    $resultsWrapper.hide();

  }

  $resultsContainer.on('click keyup', '.' + _searchItemClass, function (e) {
    if (e.type !== 'click' && e.which !== 13 && e.which !== 32) return;

    e.stopPropagation();
    e.preventDefault();

    var href = $(this).data('link') || null;

    if (href) window.location.href = href;

    _resetSearch();
  });

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
        _isSearchActive = false;
        return;
      }

      var results = index.search(query);
      var resultsHTML = '';
      var resultLength = results.length || null;

      if (query.length >= _VALID_QUERY_LENGTH && results.length === 0) {
        resultLength = 0;
      }

      if (resultLength !== null) {
        _isSearchActive = true;
        $searchContainer.addClass('search-active');
        $resultsWrapper.show();
        $searchContentBody.html('<strong><i class="fa fa-file-o"></i> ' + resultLength + '</strong> results found for <strong>' + query + '</strong>');
      } else {
        $searchContainer.removeClass('search-active');
        _isSearchActive = false;
      }

      if (results.length === 0) {
        if (_isSearchActive) {
          $searchContentBody.html('<strong>' + results.length + '</strong> results found for <strong>' + query + '</strong>');
        }
        else {
          $searchContainer.removeClass('search-active');
        }
      } else {
        for (var i = 0; i < results.length; i++) {
          var result = results[i];
          doc = documents[result.ref];
          doc.base_url = base_url;
          doc.summary = doc.text.substring(0, 200);

          resultsHTML += '<div class="' + _searchItemClass + ' animated fadeInTop" tabindex="0" role="button" data-link="' + doc.location + '">' +
            '<div class="search-body">' +
            '<h3>' + doc.title + '</h3>' +
            '<p>' + doc.summary + '</p>' +
            '</div></div>';
        }

        $results.append(resultsHTML);
        $results.highlight(query);

        setTimeout(function () { $('.' + _searchItemClass).removeClass('animated fadeInTop'); }, 500);
      }
    }

    searchInput.addEventListener('keyup', _debounce(_search, 300));

    $cancelButton.on('click', function () {
      $inputBox.val('');
      _search();
      $cancelButton.hide();
    })
  });

});