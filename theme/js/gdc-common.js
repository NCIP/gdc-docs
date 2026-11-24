$(function() {

  $(".dropdown-menu a:contains('pdf'), .bs-sidenav a:contains('pdf')").each(function () {
    var el = $(this);
    var text = el.text();
    var split = text.trim().split(' ');
    var url = split.pop();
    el.attr('href', url);
    el.attr('target', '_blank');
    el.text(split.join(' '));
    el.parent().removeClass('active');
  });

  $(".dropdown-menu a:contains('fa-'), .bs-sidenav a:contains('fa-')").each(function () {
    var el = $(this);
    var text = el.text();
    var actualContent = text.trim().split(' ').slice(1, Infinity).join(' ');
    var iconClass = text.match(/fa-\S*/);

    if (iconClass) {
      el.html(
        "<span aria-hidden='true' class='nav-icon fa " + iconClass + "'></span>"
       + actualContent
      )
    }
  });


  function ModalSearchManager(id) {
    var _self = this,
      _modalID = id,
      _modalEl,
      _modelTitleEl,
      _modelBodyTextEl,
      _searchItemClass = 'search-item';


    function _initSearch() {

      //
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
      ////////////////////////////////////////////

      //
      function _addSearchListeners() {

        var $search = $('.searchbox', _modalEl),
          submitIcon = $('.searchbox-icon', _modalEl),
          searchIconEl = submitIcon.find('.search-bttn'),
          searchCancelIconEl = submitIcon.find('.search-cancel-bttn');

        function __switchToCancelIcon() {
          searchIconEl.hide();
          searchCancelIconEl.show();
        }

        function __abortQuery(e) {
          if (e) {
            e.stopPropagation();
          }

          searchCancelIconEl.hide();
          searchIconEl.show();
          $resultsContainer.hide();
          $body.show();
          submitIcon.click();
        }

        $search.submit(function (e) {
          e.preventDefault();
          submitIcon.click();
        });

        searchCancelIconEl.click(__abortQuery);

        submitIcon.click(function () {
          if (_isSearchActive == false) {
            $inputBox.val('');
            $search.addClass('searchbox-open');
            $inputBox.focus();
            _isSearchActive = true;
          }
          else {
            $search.removeClass('searchbox-open');
            $inputBox.focusout();
            _isSearchActive = false;
          }
        });

        submitIcon.mouseup(function () {
          return false;
        });

        $search.mouseup(function () {
          return false;
        });

        $(document).mouseup(function () {
          if (_isSearchActive == true) {
            var query = $.trim($inputBox.val());

            if (query.length >= _VALID_QUERY_LENGTH) {
              __switchToCancelIcon();
            }
            else {
              __abortQuery();
            }
          }
        });
      }
      ////////////////////////////////////////////

      //
      function _getSearchTerm() {
        var sPageURL = window.location.search.substring(1);
        var sURLVariables = sPageURL.split('&');
        for (var i = 0; i < sURLVariables.length; i++) {
          var sParameterName = sURLVariables[i].split('=');
          if (sParameterName[0] == 'q') {
            return decodeURIComponent(sParameterName[1].replace(/\+/g, '%20'));
          }
        }
      }
      ////////////////////////////////////////////

      function _initSearchIndex() {
        var $results = $('.search-results', _modalEl),
          $searchContentBody = $('.search-body', _modalEl);


        $resultsContainer.delegate('.' + _searchItemClass, 'click keyup', function(e) {

          if (e.type !== 'click' && e.which !== 13 && e.which !== 32) {
            return;
          }

          e.stopPropagation();
          e.preventDefault();

          var href = $(this).data('link') || null;

          if (href) {
            window.location.href = href;
          }

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

          var documents = {},
            doc;

          for (var i = 0; i < data.docs.length; i++) {
            doc = data.docs[i];
            doc.location = base_url + doc.location;
            index.add(doc);
            documents[doc.location] = doc;
          }

          function _search() {
            var query = $inputBox.val();
            $results.empty();

            if (query.length < _VALID_QUERY_LENGTH || query === '') {
              $resultsContainer.hide();

              return;
            }

            var results = index.search(query),
                resultsHTML = '',
                resultLength = results.length || null;

            if (query.length >= _VALID_QUERY_LENGTH && results.length === 0) {
              resultLength = 0;
            }

            if (resultLength !== null) {
              $resultsContainer.show();
              $searchContentBody.html('<strong><i class="fa fa-file-o"></i> ' + resultLength  + '</strong> results found for <strong>' + query  + '</strong>' );
            }

            if (results.length > 0) {

              var baseHostURL = location.protocol + '//' + location.hostname + (location.port &&
                                                                                (location.port != 80 && location.port != 443) ? (':' + location.port) : '') +
                                '/';

              for (var i = 0; i < results.length; i++) {
                var result = results[i];
                doc = documents[result.ref];
                doc.base_url = base_url;
                doc.summary = doc.text.substring(0, 200);
                var hostURL = baseHostURL + doc.location.replace(/[\.]+\//g, '');


                resultsHTML += '<div class="' + _searchItemClass + ' animated fadeInLeft" tabindex="0" role="button" data-link="' + doc.location + '">' +
                               '<div class="doc-type-icon-container"><i class="fa fa-files-o fa-2x"></i></div>' +
                               '<div class="search-body">' +
                               '' + doc.title + '' +
                               '<p class="location-field">'  + hostURL + '&nbsp;<span class="icon-share-1"></span></p>' +
                               '<p>' + doc.summary + '</p>' +
                               '</div>' +

                               '</div>';
                //console.log(query);

              }

              $results.append(resultsHTML);
              $results.highlight(query);

              setTimeout(function() {$('.' + _searchItemClass).removeClass('animated fadeInLeft'); }, 500);
            }
            else {
              if (! _isSearchActive) {
                $body.show();
              }
              else {
                $searchContentBody.html('<strong>' + results.length  + '</strong> results found for <strong>' + query  + '</strong>' );
              }
            }
          }



          var searchInput = document.getElementById('gdc-search-query');

          var term = _getSearchTerm();

          if (term) {
            searchInput.value = term;
            search();
          }

          searchInput.addEventListener('keyup', _debounce(_search, 300));
        });

      }

      function _resetSearch() {
        $inputBox.val('');
        _self.show(false);
        $resultsContainer.hide();
      }



      function _init() {
        _addSearchListeners();
        _initSearchIndex();
      }

      var $body = $('#body'),
          $resultsContainer = $('.search-results-container', _modalEl),
          $inputBox = $('.searchbox-input', _modalEl),
          _isSearchActive = false,
          _VALID_QUERY_LENGTH = 3;


      _init();

    }


    function _init() {
      if (! _modalID) {
        console.error('Could not instantiate modal with and ID!');
        return;
      }

      _modalEl = jQuery('#' + _modalID);

      if ( _modalEl.length === 0 ) {
        console.error('Could not find modal with ID ' + _modalID +  ' !');
        return;
      }

      _modelBodyTextEl = _modalEl.find('.modal-body');

      _modalEl.on('shown.bs.modal', function() { $('#gdc-search-query').focus(); });
      _modalEl.on('hidden.bs.modal', function () { $('#gdc-search-button').focus(); });

      _initSearch();

    }

    _self.title = function(title) {

      if (arguments.length === 1) {
        _modelTitleEl.html(title);
      }

      return _modelTitleEl.text();
    };

    _self.bodyText = function(title) {

      if (arguments.length === 1) {
        _modelBodyTextEl.html(title);
      }

      return _modelBodyTextEl.text();
    };

    _self.show = function(shouldShow) {
      var toggleArg = shouldShow === false ? 'hide' : 'show';

      _modalEl.modal(toggleArg);
    };


    _init();

  }

  function init() {

    // Initialize a JS global to be used with dynamic JS Apps

    window.$gdcApp = window.$gdcApp || {};

    window.$gdcApp.config = {};

    window.$gdcApp.searchModal = new ModalSearchManager('search-modal');

    function _initScrollSpy() {

      if ($('#body').data('current-page') === 'Home') return;

      var scrollSpyTarget = '.bs-sidebar',
        scrollBody = $('html, body');



      // Stripe tables
      $('table').addClass('table table-striped table-hover');

      // Enable side ToC
      $('body').scrollspy({
        target: scrollSpyTarget,
        offset: 0
      });

      var sideBar = $('.bs-sidebar');

      $(scrollSpyTarget + ' a[href^=\'#\']').on('click', function(e) {

        // prevent default anchor click behavior
        e.preventDefault();

        // store hash
        var hash = this.hash,
          scrollTargetEl = $(hash);

          // animate
          scrollBody.animate({
            scrollTop: scrollTargetEl.offset().top + 7 /* plus some delta so scrollSpy highlight gets triggered */
          }, 300,
            function(){

            var targetEl = scrollTargetEl.find('a'),
              classes = 'animated-long focusText';

            targetEl.addClass(classes);

            setTimeout(function() { targetEl.removeClass(classes); }, 1100);
            // when done, add hash to url
            // (default click behaviour)
            window.location.hash = hash;

        });

      });

      var mainContainer = $('.main-container'),
          selectedNavRegion = sideBar.find('.main'),
          _totalAnchorHeight = 0,
          _anchorOffsetMap = [];

      sideBar.find('a').each(function() {
        var anchor = $(this),
          anchorHeight = anchor.outerHeight();

        _anchorOffsetMap.push({offset: _totalAnchorHeight, height: anchorHeight});

        _totalAnchorHeight += anchorHeight;
      });

      // TODO: Could improve the search complexity O(n) given that the list is sorted by offset
      var findAnchorForOffset = function(offset) {

        if (! _anchorOffsetMap.length) {
          return 0;
        }
        else if (_anchorOffsetMap.length === 1) {
          return _anchorOffsetMap[0];
        }

        var i = 0;

        for (; i < _anchorOffsetMap.length; i++) {
          if (_anchorOffsetMap[i].offset > offset) {
            break;
          }
        }

        return _anchorOffsetMap[i - 1];
      };

      $(window).scroll(function() {
        var scollableDistance = Math.max(0, mainContainer.outerHeight() + mainContainer.offset().top +
                                            $('#docs-footer').outerHeight() - $(window).outerHeight());
        var percentPageScrolled = Math.min(1.0, $(window).scrollTop() / scollableDistance);

        var proposedOffset = selectedNavRegion.outerHeight() * percentPageScrolled,
          anchorOffset = findAnchorForOffset(proposedOffset);

        sideBar.stop().animate({
          scrollTop: Math.min(anchorOffset.offset - anchorOffset.height * 2, proposedOffset)
        }, 200);

      });

      sideBar.scroll(function(e) { e.stopPropagation(); });

      // Prevent disabled links from causing a page reload
      $('li.disabled a').click(function (e) {
        e.preventDefault();
      });
    }



    function _initHeaders(confineToContainerID) {

      var $body = $(confineToContainerID);

      // Add header links
      $(":header", $body).each(function (i, header) {
        var $header = $(header);

        if ($header.hasClass('no-auto-render')) {
          return;
        }

        var id = $header.attr('id');
        var icon = '&nbsp;<i class="fa fa-share-alt"></i>';

        if (id) {
          var title = $header.text();
          $header.text("");
          //$header.prepend($("<a/>").addClass("header-link").attr("href", "#" + id)); //.html(icon));
          $header.append($("<a/>")
            .addClass("header-text-link")
            .attr("href", "#" + id)
            .attr("title", "Click on this header and copy URL to link to this section.")
            .append(title)
            .append(icon));
        }
      });


      var mainHeader = $('h1', $body);

      if (mainHeader.hasClass('no-auto-render')) {
        return;
      }

      mainHeader.prepend('<span class="header-badge"><span class="fa fa-book" aria-hidden="true"></span></span>');

    }

    function _initLinks(confineToContainerID) {
      // Ensure all anchor links which appear to be external have the appropriate target and
      // icon.
      $(confineToContainerID + ' a[href^="http"]').each(function() {
        var anchor = $(this);

        anchor.attr('target', '_blank');
        anchor.addClass('external-link');

        // If the link already has an font awesome icon associated skip it
        if (! anchor.find('*[class*="fa-"]').length) {
          anchor.prepend('&nbsp;<i class="fa fa-external-link"></i>&nbsp;');
        }

      });


      // Ensure all mail links have the appropriate icon.
      $('a[href^="mailto"]').prepend('&nbsp;<i class="fa fa-envelope"></i>&nbsp;');
    }

    function _initMenuNavBar(container, subContainer) {
      var menuBar = $('.menu-bar'),
        container = $(container);


      container.find(subContainer)
        .hover(function () {
            var dropdownItem = $(this);
            menuBar.css({width: dropdownItem.width(), left: dropdownItem.position().left});
          },
          function () {
            menuBar.css({width: 0});
          });

      $('.navbar.navbar-default').autoHidingNavbar({hideOffset: _hideMenuOffset});
    }

    function _ensureMaxHeight(resizingEl) {
      var windowEl = $(window),
          footer = $('#docs-footer'),
          footerHeight = footer.outerHeight(),
          footerOffsetTop = footer.offset().top,
          OFFSET = -(130 + footerHeight),
          resizeElHeight = 0;


      function _recalcMax() {
        var resizeHeight = resizingEl.outerHeight(),
            distance = Math.round(footerOffsetTop - resizeHeight - windowEl.scrollTop() - footerHeight /* <-- offset */),
            windowHeight =  windowEl.height();

        if (distance > 0) {
          resizingEl.css({overflow: 'auto', maxHeight: (windowHeight - Math.round(1.5 * _hideMenuOffset)) + 'px'}, 'fast');
        }
        else {
          resizeElHeight = windowHeight + OFFSET;
          resizingEl.css({overflow: 'auto', maxHeight: resizeElHeight + 'px'});
        }

      }

      var _prevScrollOffset = 0,
          _document = $(document);


      function _onScroll() {
        var currentScrollOffset = windowEl.scrollTop(),
            delta = currentScrollOffset - _prevScrollOffset,
            absDelta = Math.abs(delta),
            scrollPosition = windowEl.height() + currentScrollOffset,
            scrollHeight = _document.height(),
            shouldResetNav = (scrollPosition + _hideMenuOffset) >=  scrollHeight;

        if (delta > 0 && absDelta >= _hideMenuOffset && ! shouldResetNav) {
          resizingEl.css({top: 10});
          _prevScrollOffset = currentScrollOffset;
        }
        else if (delta < 0 && absDelta >= _hideMenuOffset || shouldResetNav) {
          resizingEl.css({top: ''});
          _prevScrollOffset = currentScrollOffset;
        }

        _recalcMax();
      }


      windowEl.scroll(_onScroll);
      windowEl.resize(_recalcMax)
    }



    function _calcMainContentWidth() {
      if ($('.full-width-content').length || $('.toc-container').length === 0) {
        $('.main-container').addClass('col-md-12').removeClass('col-md-9').css({borderLeft: 'none'});
      }
    }

    function _initAlerts() {
      $('.alert').each(function() {
        var alertContainer = $(this),
            supportedAlertTypes = ['alert-info'];

        for (var i = 0; i < supportedAlertTypes.length; i++) {
          var alertTypeClass = supportedAlertTypes[i];

          if (!  alertContainer.hasClass(alertTypeClass)) {
            continue;
          }

          var iconClass = null,
              alertText = '';

          switch (alertTypeClass) {
            case 'alert-info':
              iconClass = 'icon-pencil';
            break;
            case 'alert-warning':
              iconClass = 'icon-attention-1';
              break;
           default:
             break;
          }

          if( ! iconClass ) {
            continue;
          }


          alertText = alertContainer.text();
          alertContainer.empty();
          alertContainer.append('<div class="alert-indicator-icon">' +
                                '<i class="' + iconClass + '"></i>' +
                                '</div>' +
                                '<div class="alert-indicator-text">' +
                                alertText + '</div>');


        }

      });
    }

    function _initScrollUpIndicator() {
      $.scrollUp({
        scrollName: 'scroll-up-indicator',
        scrollDistance: 300,
        scrollFrom: 'top',
        scrollSpeed: 300,
        easingType: 'swing',
        animation: 'fade',
        animationSpeed: 200,
        scrollText: 'Scroll to top',
        scrollTitle: 'Scroll to the top of this page.',
        scrollImg: {
          active: true
        },
        activeOverlay: false,
        zIndex: 100000
      });

      $('#scroll-up-indicator').html('<span style="display: none">Scroll to the top of this page.</span>');
    }

    var _bsSidebar = $('.bs-sidebar');

    if (_bsSidebar.length) {
      _ensureMaxHeight(_bsSidebar);
    }

    var BODY_ID = '#body';
    var _hideMenuOffset = 64;

    _initScrollSpy();
    _initMenuNavBar('.navbar-nav', '> li');
    _initHeaders(BODY_ID);
    _initLinks(BODY_ID);
    _calcMainContentWidth();
    _initAlerts();
    _initScrollUpIndicator();


    // Hightlight code
    hljs.initHighlightingOnLoad();
    /*
    var particleContainer = $('.spinParticle');
    particleContainer.addClass('endLoad');

    setTimeout(function () {
      particleContainer.hide();
    }, 1000);

    var _handleFontTransition = function () {
      var bodyEl = $('.main-container');

      if (bodyEl.hasClass('loading-content')) {
        bodyEl.removeClass('fadeInBlurIntro loading-content').addClass('fadeInBlur');
      }

      $('.spinParticle').fadeOut('fast');
    };

    setTimeout(function() {
      fontSpy('FontAwesome', {
        glyphs: '\ue800\ue8019\ue81c\ue843',
        success: _handleFontTransition,
        failure: _handleFontTransition
      });

    }, 0);*/




  }

  /////////////////////////////////////////
  init();

});
