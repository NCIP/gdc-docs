(function (fetch, d3, Remarkable) {


  var getNotifications = function() {
    var checkStatus = function(response) {
      if (response.status >= 200 && response.status < 300) {
        return response;
      } else {
        var error = new Error(response.statusText);
        error.response = response;
        throw error;
      }
    };

    return fetch('https://api.gdc.cancer.gov/notifications')
      .then(checkStatus)
      .then(function(response) {
        return response.json();
      }).then(function(json) {
        return json.data || [];
      });
  };

  var readDismissedCookie = function() {
    var cookieValue = document.cookie.replace(/(?:(?:^|.*;\s*)gdc-dismissed-notifications\s*\=\s*([^;]*).*$)|^.*$/, "$1");
    var dismissedIds = cookieValue ? JSON.parse(cookieValue) : [];
    return dismissedIds;
  };

  var addDismissedIdToCookie = function(notificationID) {
    var dismissedIds = readDismissedCookie();
    if (dismissedIds.indexOf(notificationID) === -1) {
      dismissedIds = dismissedIds.concat(notificationID);
      var cookieString = 'gdc-dismissed-notifications='
        .concat(JSON.stringify(dismissedIds))
        .concat(';path=/');
      document.cookie=cookieString;
    }
  };

  var renderNotifications = function(notifications) {
    var container = d3.select('#notifications');

    var markdownParser = new Remarkable();
    var bannerHTML = '<span class="fa icon"> \
      </span>\
      <span class="header-message"> \
      </span> \
      <span class="header-banner-dismiss"> \
          Dismiss <span class="fa fa-close" aria-hidden="true"></span>\
      </span>';

    var BANNER_HEIGHT = 40;
    d3.select('#docs-container').style('padding-top', notifications.length * BANNER_HEIGHT + 'px');
    d3.select('.bs-sidebar').style('margin-top', notifications.length * BANNER_HEIGHT + 'px');
    notifications.forEach(function(n) {
      var banner = container.append('div')
      .classed('header-banner', true)
      .classed('enter', true)
      .html(bannerHTML);
      setTimeout(function() { banner.classed('enter', false) });
      banner.classed('warning', n.level === 'WARNING')
            .classed('error', n.level === 'ERROR');

      banner.select('.fa-circle-o')
      .classed('hidden', n.level === 'ERROR');

      banner.select('.icon')
      .classed('fa-question', n.level === 'INFO')
      .classed('fa-exclamation', n.level === 'WARNING')
      .classed('fa-exclamation-triangle', n.level === 'ERROR');

      banner.select('.header-message')
            .html(markdownParser.render(n.message));

      banner.select('.header-banner-dismiss')
            .classed('hidden', !n.dismissible)
            .on('click', function() {
              banner.classed('dismissed', n.dismissible);
              var paddingTop = d3.select('#docs-container').style('padding-top').replace(/px/, '');
              d3.select('#docs-container').style('padding-top', parseInt(paddingTop, 10) - BANNER_HEIGHT + 'px');
              d3.select('.bs-sidebar').style('margin-top', parseInt(paddingTop, 10) - BANNER_HEIGHT + 'px');
              addDismissedIdToCookie(n.id);
            });

    });

  }

  getNotifications().then(function(data) {
    var alreadyDismissedIds = readDismissedCookie();
    var docNotifications = data
      .filter(function(d) { return (d.components || []).indexOf('DOCUMENTATION') !== -1; })
      .filter(function(d) { return alreadyDismissedIds.indexOf(d.id) === -1; } )
      .reduce((sorted, d) => d.dismissible ? sorted.concat(d) : [d].concat(sorted), []);

    if (docNotifications.length) {
      renderNotifications(docNotifications);
    }
  });

})(window.fetch, d3, Remarkable);
