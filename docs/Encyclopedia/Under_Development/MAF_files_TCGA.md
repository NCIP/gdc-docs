# Mutation Annotation Format (MAF) Specification - v2.4 #
Draft
## Description ##
##Current version changes##
This current revision is version 2.4 of the Mutation Annotation Format (MAF) specification.

The following items in the specification were added or modified in version 2.4 from version 2.3:

Header for MAF file is "#version 2.4"
Allowed "Validation_Status" values are now "Untested, Inconclusive, Valid, and Invalid."
"Mutation_Status" values are now aligned with the VCF VLS field and are "None, Germline, Somatic, LOH, Post-transcriptional modification, Unknown" (addition of the enumerated value of "Post-transcriptional modification")
The VCF VLS value of Wild Type will be replaced by None
The value in the "Validation_Status" field determines what values are allowed in the "Mutation_Status" field
The "Validation_Status" field will no longer accept NULL values
The values allowed in the "Sequence_Source" column have been changed to be a subset of the SRA 1.5 library_strategy values.
All Version 2.4 MAF files will be required to contain UUID columns (Tumor_Sample_UUID and Matched_Norm_Sample_UUID)
For the "Variant_Classification" field, the values of De_novo_Start_InFrame and De_novo_Start_OutOfFrame are no longer allowed.
The values for the"dbSNP_Val_Status" will be enforced. "none" will no longer be allowed.
"Somatic" is the only acceptable value for 'Mutation_Status' for a somatic MAF (named .somatic.maf). Protected MAF (named .protected.maf) has no such restriction and can contain Somatic, Germline, Unknown, LOH, Post-transcriptional modification, None for Mutation_Status.
For a somatic MAF, following rules should be satisfied:
SOMATIC = A AND (B OR C OR D)
A: Mutation_Status == "Somatic"
B: Validation_Status == "Valid"
C. Verification_Status == "Verified"
D. Variant_Classification is not {Intron, 5'UTR, 3'UTR, 5'Flank, 3'Flank, IGR}, which implies that Variant_Classification can only be \{Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region}.
If Validation_Status  == "Valid" then Tumor_Validation_Allele1, Tumor_Validation_Allele2, Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2 cannot  be null
MAF 2.3 Validation rule 8a was superceded by:
If Validation_Status == "Invalid" then Tumor_Validation_Allele1, Tumor_Validation_Allele2, Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2 cannot be null AND Tumor_Validation_Allelle1 == Match_Norm_Validation_Allele1 AND Tumor_Validation_Allelle2 == Match_Norm_Validation_Allele2  (Added as a replacement for 8a as a result of breakdown)


## About MAF specifications ##
Mutation annotation files should be transferred to the DCC. Those files should be formatted using the mutation annotation format (MAF) that is described below. File naming convention is also below.

Following categories of somatic mutations are reported in MAF files:

* Missense and nonsense
* Splice site, defined as SNP within 2 bp of the splice junction
* Silent mutations
* Indels that overlap the coding region or splice site of a gene or the targeted region of a genetic element of interest.
* Frameshift mutations
* Mutations in regulatory regions




Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki       var contextPath = '';     .ia-fixed-sidebar, .ia-splitter-left {width: 285px;}.theme-default .ia-splitter #main {margin-left: 285px;}.ia-fixed-sidebar {visibility: hidden;}             <fieldset class="i18n hidden"><input type="hidden" name="calendar3.month.long.july" value="July"><input type="hidden" name="calendar3.day.short.wednesday" value="Wed"><input type="hidden" name="calendar3.day.short.thursday" value="Thu"><input type="hidden" name="calendar3.month.short.march" value="Mar"><input type="hidden" name="calendar3.month.long.april" value="April"><input type="hidden" name="calendar3.month.long.october" value="October"><input type="hidden" name="calendar3.month.long.august" value="August"><input type="hidden" name="calendar3.month.short.july" value="Jul"><input type="hidden" name="calendar3.month.short.may" value="May"><input type="hidden" name="calendar3.month.short.november" value="Nov"><input type="hidden" name="calendar3.day.long.friday" value="Friday"><input type="hidden" name="calendar3.day.long.sunday" value="Sunday"><input type="hidden" name="calendar3.day.long.saturday" value="Saturday"><input type="hidden" name="calendar3.month.short.april" value="Apr"><input type="hidden" name="calendar3.day.long.wednesday" value="Wednesday"><input type="hidden" name="calendar3.month.long.december" value="December"><input type="hidden" name="calendar3.month.short.october" value="Oct"><input type="hidden" name="calendar3.day.long.monday" value="Monday"><input type="hidden" name="calendar3.month.short.june" value="Jun"><input type="hidden" name="calendar3.day.short.monday" value="Mon"><input type="hidden" name="calendar3.day.short.tuesday" value="Tue"><input type="hidden" name="calendar3.day.short.saturday" value="Sat"><input type="hidden" name="calendar3.month.long.march" value="March"><input type="hidden" name="calendar3.month.long.june" value="June"><input type="hidden" name="calendar3.month.short.february" value="Feb"><input type="hidden" name="calendar3.month.short.august" value="Aug"><input type="hidden" name="calendar3.month.short.december" value="Dec"><input type="hidden" name="calendar3.day.short.sunday" value="Sun"><input type="hidden" name="calendar3.month.long.february" value="February"><input type="hidden" name="calendar3.day.long.tuesday" value="Tuesday"><input type="hidden" name="calendar3.month.long.may" value="May"><input type="hidden" name="calendar3.month.long.september" value="September"><input type="hidden" name="calendar3.month.long.november" value="November"><input type="hidden" name="calendar3.month.short.january" value="Jan"><input type="hidden" name="calendar3.month.short.september" value="Sep"><input type="hidden" name="calendar3.day.long.thursday" value="Thursday"><input type="hidden" name="calendar3.month.long.january" value="January"><input type="hidden" name="calendar3.day.short.friday" value="Fri"></fieldset>      <div class="gliffy-webpanel-footer"><span></span></div>                                        window.WRM=window.WRM||{};window.WRM.\_unparsedData=window.WRM.\_unparsedData||{};window.WRM.\_unparsedErrors=window.WRM.\_unparsedErrors||{}; WRM.\_unparsedData\["com.atlassian.plugins.atlassian-plugins-webresource-plugin:context-path.context-path"\]="\\u0022\\u0022"; WRM.\_unparsedData\["com.atlassian.confluence.plugins.confluence-hipchat-integration-plugin:discovery-javascript-data.link-active"\]="{\\u0022linkActive\\u0022:false,\\u0022conditionsMet\\u0022:false,\\u0022admin\\u0022:false}"; WRM.\_unparsedData\["com.atlassian.applinks.applinks-plugin:applinks-common-exported.applinks-help-paths"\]="{\\u0022entries\\u0022:{\\u0022applinks.docs.root\\u0022:\\u0022https://confluence.atlassian.com/display/APPLINKS-052/\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.sslunmatched\\u0022:\\u0022SSL+and+application+link+troubleshooting+guide\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.oauthsignatureinvalid\\u0022:\\u0022OAuth+troubleshooting+guide\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.oauthtimestamprefused\\u0022:\\u0022OAuth+troubleshooting+guide\\u0022,\\u0022applinks.docs.delete.entity.link\\u0022:\\u0022Create+links+between+projects\\u0022,\\u0022applinks.docs.adding.application.link\\u0022:\\u0022Link+Atlassian+applications+to+work+together\\u0022,\\u0022applinks.docs.administration.guide\\u0022:\\u0022Application+Links+Documentation\\u0022,\\u0022applinks.docs.oauth.security\\u0022:\\u0022OAuth+security+for+application+links\\u0022,\\u0022applinks.docs.troubleshoot.application.links\\u0022:\\u0022Troubleshoot+application+links\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.unknownerror\\u0022:\\u0022Network+and+connectivity+troubleshooting+guide\\u0022,\\u0022applinks.docs.configuring.auth.trusted.apps\\u0022:\\u0022Configuring+Trusted+Applications+authentication+for+an+application+link\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.authlevelunsupported\\u0022:\\u0022OAuth+troubleshooting+guide\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.ssluntrusted\\u0022:\\u0022SSL+and+application+link+troubleshooting+guide\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.unknownhost\\u0022:\\u0022Network+and+connectivity+troubleshooting+guide\\u0022,\\u0022applinks.docs.delete.application.link\\u0022:\\u0022Link+Atlassian+applications+to+work+together\\u0022,\\u0022applinks.docs.link.applications\\u0022:\\u0022Link+Atlassian+applications+to+work+together\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.oauthproblem\\u0022:\\u0022OAuth+troubleshooting+guide\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.migration\\u0022:\\u0022Update+application+links+to+use+OAuth\\u0022,\\u0022applinks.docs.relocate.application.link\\u0022:\\u0022Link+Atlassian+applications+to+work+together\\u0022,\\u0022applinks.docs.administering.entity.links\\u0022:\\u0022Create+links+between+projects\\u0022,\\u0022applinks.docs.upgrade.application.link\\u0022:\\u0022OAuth+security+for+application+links\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.connectionrefused\\u0022:\\u0022Network+and+connectivity+troubleshooting+guide\\u0022,\\u0022applinks.docs.configuring.auth.oauth\\u0022:\\u0022OAuth+security+for+application+links\\u0022,\\u0022applinks.docs.insufficient.remote.permission\\u0022:\\u0022OAuth+security+for+application+links\\u0022,\\u0022applinks.docs.configuring.application.link.auth\\u0022:\\u0022OAuth+security+for+application+links\\u0022,\\u0022applinks.docs.diagnostics\\u0022:\\u0022Application+links+diagnostics\\u0022,\\u0022applinks.docs.configured.authentication.types\\u0022:\\u0022OAuth+security+for+application+links\\u0022,\\u0022applinks.docs.adding.entity.link\\u0022:\\u0022Create+links+between+projects\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.unexpectedresponse\\u0022:\\u0022Network+and+connectivity+troubleshooting+guide\\u0022,\\u0022applinks.docs.configuring.auth.basic\\u0022:\\u0022Configuring+Basic+HTTP+Authentication+for+an+Application+Link\\u0022,\\u0022applinks.docs.diagnostics.troubleshoot.authlevelmismatch\\u0022:\\u0022OAuth+troubleshooting+guide\\u0022}}"; WRM.\_unparsedData\["com.atlassian.applinks.applinks-plugin:applinks-common-exported.applinks-types"\]="{\\u0022crowd\\u0022:\\u0022Crowd\\u0022,\\u0022confluence\\u0022:\\u0022Confluence\\u0022,\\u0022fecru\\u0022:\\u0022FishEye / Crucible\\u0022,\\u0022stash\\u0022:\\u0022Bitbucket Server\\u0022,\\u0022jira\\u0022:\\u0022JIRA\\u0022,\\u0022refapp\\u0022:\\u0022Reference Application\\u0022,\\u0022bamboo\\u0022:\\u0022Bamboo\\u0022,\\u0022generic\\u0022:\\u0022Generic Application\\u0022}"; WRM.\_unparsedData\["com.atlassian.applinks.applinks-plugin:applinks-common-exported.entity-types"\]="{\\u0022singular\\u0022:{\\u0022refapp.charlie\\u0022:\\u0022Charlie\\u0022,\\u0022fecru.project\\u0022:\\u0022Crucible Project\\u0022,\\u0022fecru.repository\\u0022:\\u0022FishEye Repository\\u0022,\\u0022stash.project\\u0022:\\u0022Bitbucket Server Project\\u0022,\\u0022generic.entity\\u0022:\\u0022Generic Project\\u0022,\\u0022confluence.space\\u0022:\\u0022Confluence Space\\u0022,\\u0022bamboo.project\\u0022:\\u0022Bamboo Project\\u0022,\\u0022jira.project\\u0022:\\u0022JIRA Project\\u0022},\\u0022plural\\u0022:{\\u0022refapp.charlie\\u0022:\\u0022Charlies\\u0022,\\u0022fecru.project\\u0022:\\u0022Crucible Projects\\u0022,\\u0022fecru.repository\\u0022:\\u0022FishEye Repositories\\u0022,\\u0022stash.project\\u0022:\\u0022Bitbucket Server Projects\\u0022,\\u0022generic.entity\\u0022:\\u0022Generic Projects\\u0022,\\u0022confluence.space\\u0022:\\u0022Confluence Spaces\\u0022,\\u0022bamboo.project\\u0022:\\u0022Bamboo Projects\\u0022,\\u0022jira.project\\u0022:\\u0022JIRA Projects\\u0022}}"; WRM.\_unparsedData\["com.atlassian.applinks.applinks-plugin:applinks-common-exported.authentication-types"\]="{\\u0022com.atlassian.applinks.api.auth.types.BasicAuthenticationProvider\\u0022:\\u0022Basic Access\\u0022,\\u0022com.atlassian.applinks.api.auth.types.TrustedAppsAuthenticationProvider\\u0022:\\u0022Trusted Applications\\u0022,\\u0022com.atlassian.applinks.api.auth.types.CorsAuthenticationProvider\\u0022:\\u0022CORS\\u0022,\\u0022com.atlassian.applinks.api.auth.types.OAuthAuthenticationProvider\\u0022:\\u0022OAuth\\u0022,\\u0022com.atlassian.applinks.api.auth.types.TwoLeggedOAuthAuthenticationProvider\\u0022:\\u0022OAuth\\u0022,\\u0022com.atlassian.applinks.api.auth.types.TwoLeggedOAuthWithImpersonationAuthenticationProvider\\u0022:\\u0022OAuth\\u0022}"; WRM.\_unparsedData\["com.atlassian.confluence.plugins.confluence-feature-discovery-plugin:confluence-feature-discovery-plugin-resources.test-mode"\]="false"; WRM.\_unparsedData\["com.atlassian.confluence.plugins.confluence-license-banner:confluence-license-banner-resources.license-details"\]="{\\u0022daysBeforeLicenseExpiry\\u0022:0,\\u0022daysBeforeMaintenanceExpiry\\u0022:0,\\u0022showLicenseExpiryBanner\\u0022:false,\\u0022showMaintenanceExpiryBanner\\u0022:false,\\u0022renewUrl\\u0022:null,\\u0022salesEmail\\u0022:null}"; if(window.WRM.\_dataArrived)window.WRM.\_dataArrived();                    function recordOutboundLink(link, category, action) { try { var myTracker=\_gat.\_getTrackerByName(); \_gaq.push(\['myTracker.\_trackEvent', ' + category + ', ' + action + '\]); setTimeout('document.location = "' + link.href + '"', 100) }catch(err){} }   // Standard Google Universal Analytics code (function(i,s,o,g,r,a,m){i\['GoogleAnalyticsObject'\]=r;i\[r\]=i\[r\]||function(){ (i\[r\].q=i\[r\].q||\[\]).push(arguments)},i\[r\].l=1*new Date();a=s.createElement(o), m=s.getElementsByTagName(o)\[0\];a.async=1;a.src=g;m.parentNode.insertBefore(a,m) })(window,document,'script','//www.google-analytics.com/analytics.js','ga'); AJS.toInit(function(){ ga('create', 'UA-79015079-1', 'auto'); // Add a page-level custom variable to record the space-key if (typeof AJS.params.spaceKey === 'string') { ga('set', 'dimension1', AJS.params.spaceKey); // Set a \`spaceKey\` dimension at page level } ga('send', 'pageview'); });    // Standard Google Universal Analytics code (function(i,s,o,g,r,a,m){i\['GoogleAnalyticsObject'\]=r;i\[r\]=i\[r\]||function(){ (i\[r\].q=i\[r\].q||\[\]).push(arguments)},i\[r\].l=1*new Date();a=s.createElement(o), m=s.getElementsByTagName(o)\[0\];a.async=1;a.src=g;m.parentNode.insertBefore(a,m) })(window,document,'script','//www.google-analytics.com/analytics.js','ga'); AJS.toInit(function(){ ga('create', 'UA-114974387-1', 'auto'); // Add a page-level custom variable to record the space-key if (typeof AJS.params.spaceKey === 'string') { ga('set', 'dimension1', AJS.params.spaceKey); // Set a \`spaceKey\` dimension at page level } ga('send', 'pageview'); });   .expandable-text-data .expandable-text-data-content {max-height: undefined}                                         .tableFloatingHeader{display:none !important;}.tableFloatingHeaderOriginal{position:static !important;}.tableFloatingHeader{display:none !important;}.tableFloatingHeaderOriginal{position:static !important;}

[![Skip Navigation](./Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki_files/skipnav.gif "Skip Navigation")](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#maincontent)



*   [Skip to content](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#title-heading)
*   [Skip to breadcrumbs](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#breadcrumbs)
*   [Skip to header menu](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#header-menu-bar)
*   [Skip to action menu](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#navigation)
*   [Skip to quick search](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#quick-search-query)

#preview #previewArea iframe { overflow: visible; }

[![NIH | National Cancer Institute | NCI Wiki](./Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki_files/NCI_Wiki_54px_Logo_COLOR_cropped.png "NIH | National Cancer Institute | NCI Wiki")](http://www.cancer.gov/)

 

[![New Account](./Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki_files/new_account.png "New Account")](http://wikiutils.nci.nih.gov/wiki_signup)

[![Help Tips](./Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki_files/helpful_tips.png "Help Tips")](https://wiki.nci.nih.gov/x/-wxy)

[Linked Applications](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#app-switcher)

*   [NCI Tracker](https://tracker.nci.nih.gov/secure/MyJiraHome.jspa "https://tracker.nci.nih.gov/secure/MyJiraHome.jspa")
*   [Collaborate](https://collaborate.nci.nih.gov/ "https://collaborate.nci.nih.gov/")
*   [National Cancer Institute - Confluence Wiki](https://wiki.nci.nih.gov/ "https://wiki.nci.nih.gov/")

[National Cancer Institute - Confluence Wiki](https://wiki.nci.nih.gov/)
========================================================================

*   [Spaces](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#space-menu-link-content "Spaces")

*   [People](https://wiki.nci.nih.gov/browsepeople.action "People")
*   [Calendars](https://wiki.nci.nih.gov/calendar/mycalendar.action "Calendars")
*   [Create](https://wiki.nci.nih.gov/pages/createpage.action?spaceKey=TCGA&fromPageId=26187384&src=quick-create "Create blank page (c)") [Create](https://wiki.nci.nih.gov/pages/createpage.action?spaceKey=TCGA&fromPageId=26187384 "Create from template")

*   Quick Search

*   [Help](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "Help")

   *   [Online Help](https://docs.atlassian.com/confluence/docs-510/ "Visit the Confluence documentation home")
   *   [Keyboard Shortcuts](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "View available keyboard shortcuts (?)")
   *   [Feed Builder](https://wiki.nci.nih.gov/dashboard/configurerssfeed.action "Create your custom RSS feed.")
   *   [Available Gadgets](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "Browse gadgets provided by Confluence")
   *   [About Confluence](https://wiki.nci.nih.gov/aboutconfluencepage.action "Get more information about Confluence")


*   [

   0](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "Open Notifications (g , n)")
*   [

   ![](./Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki_files/default.png)

   ](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "Miller, Ian (Federation/Google)")

   *   [Recently Viewed](https://wiki.nci.nih.gov/users/viewuserhistory.action " (g , r)")

   *   [Profile](https://wiki.nci.nih.gov/users/viewmyprofile.action)
   *   [Tasks](https://wiki.nci.nih.gov/plugins/inlinetasks/mytasks.action)
   *   [Saved for later](https://wiki.nci.nih.gov/users/viewmyfavourites.action)
   *   [Watches](https://wiki.nci.nih.gov/users/viewnotifications.action)
   *   [Drafts](https://wiki.nci.nih.gov/users/viewmydrafts.action)
   *   [Network](https://wiki.nci.nih.gov/users/viewfollow.action?username=milleri)
   *   [Settings](https://wiki.nci.nih.gov/users/viewmysettings.action)
   *   [Atlassian Marketplace](https://wiki.nci.nih.gov/plugins/servlet/upm/requests?source=header_user)

   *   [Log Out](https://wiki.nci.nih.gov/logout.action)




[![TCGA](./Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki_files/default-space-logo-256.png)](https://wiki.nci.nih.gov/display/TCGA/TCGA+Wiki+Home?src=sidebar "TCGA")

[TCGA](https://wiki.nci.nih.gov/display/TCGA/TCGA+Wiki+Home?src=sidebar "TCGA")

*   [Pages](https://wiki.nci.nih.gov/collector/pages.action?key=TCGA&src=sidebar)

##### Page tree

*   [TCGA Home](https://wiki.nci.nih.gov/display/TCGA/TCGA+Home?src=contextnavpagetreemode)

*   [About TCGA](https://wiki.nci.nih.gov/display/TCGA/About+TCGA?src=contextnavpagetreemode)

*   [About TCGA Data](https://wiki.nci.nih.gov/display/TCGA/About+TCGA+Data?src=contextnavpagetreemode)

*   [Access Tiers](https://wiki.nci.nih.gov/display/TCGA/Access+Tiers?src=contextnavpagetreemode)

*   [Data Levels and Data Types](https://wiki.nci.nih.gov/display/TCGA/Data+Levels+and+Data+Types?src=contextnavpagetreemode)

*   [Platform Design](https://wiki.nci.nih.gov/display/TCGA/Platform+Design?src=contextnavpagetreemode)

*   [TCGA Encyclopedia](https://wiki.nci.nih.gov/display/TCGA/TCGA+Encyclopedia?src=contextnavpagetreemode)

*   [Software Release Documents](https://wiki.nci.nih.gov/display/TCGA/Software+Release+Documents?src=contextnavpagetreemode)

*   [Applications](https://wiki.nci.nih.gov/display/TCGA/Applications?src=contextnavpagetreemode)

*   [Data Reports and Dashboards](https://wiki.nci.nih.gov/display/TCGA/Data+Reports+and+Dashboards?src=contextnavpagetreemode)

*   [Web Services](https://wiki.nci.nih.gov/display/TCGA/Web+Services?src=contextnavpagetreemode)

*   [TCGA User Documentation](https://wiki.nci.nih.gov/display/TCGA/TCGA+User+Documentation?src=contextnavpagetreemode)

   *   [Publication Guidelines](https://wiki.nci.nih.gov/display/TCGA/Publication+Guidelines?src=contextnavpagetreemode)

   *   [TCGA Specifications](https://wiki.nci.nih.gov/display/TCGA/TCGA+Specifications?src=contextnavpagetreemode)

       *   [File Format Specifications](https://wiki.nci.nih.gov/display/TCGA/File+Format+Specifications?src=contextnavpagetreemode)

           *   [Biospecimen and Clinical XSD Files Specification](https://wiki.nci.nih.gov/display/TCGA/Biospecimen+and+Clinical+XSD+Files+Specification?src=contextnavpagetreemode)

           *   [GSC MAGE-TAB Specification](https://wiki.nci.nih.gov/display/TCGA/GSC+MAGE-TAB+Specification?src=contextnavpagetreemode)

           *   [Mutation Annotation Format (MAF) Specification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification?src=contextnavpagetreemode)

               *   [Mutation Annotation Format (MAF) Specification - v2.4](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4?src=contextnavpagetreemode)

               *   [Mutation Annotation Format (MAF) Specification - v2.3](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.3?src=contextnavpagetreemode)

               *   [Mutation Annotation Format (MAF) Specification - v2.2](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.2?src=contextnavpagetreemode)

               *   [Mutation Annotation Format (MAF) Specification - v2.1](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.1?src=contextnavpagetreemode)

               *   [Mutation Annotation Format (MAF) Specification - v2.0](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.0?src=contextnavpagetreemode)

               *   [Mutation Annotation Format (MAF) Specification - v1.0](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v1.0?src=contextnavpagetreemode)


           *   [Protein Array Data Format Specification](https://wiki.nci.nih.gov/display/TCGA/Protein+Array+Data+Format+Specification?src=contextnavpagetreemode)

           *   [RNASeq Data Format Specification](https://wiki.nci.nih.gov/display/TCGA/RNASeq+Data+Format+Specification?src=contextnavpagetreemode)

           *   [Wiggle Format Specification](https://wiki.nci.nih.gov/display/TCGA/Wiggle+Format+Specification?src=contextnavpagetreemode)

           *   [TCGA Variant Call Format (VCF) Specification](https://wiki.nci.nih.gov/display/TCGA/TCGA+Variant+Call+Format+%28VCF%29+Specification?src=contextnavpagetreemode)



   *   [TCGA User's Guides](https://wiki.nci.nih.gov/display/TCGA/TCGA+User%27s+Guides?src=contextnavpagetreemode)


*   [External Analytical Tools](https://wiki.nci.nih.gov/display/TCGA/External+Analytical+Tools?src=contextnavpagetreemode)

*   [Frequently Asked Questions (FAQ)](https://wiki.nci.nih.gov/pages/viewpage.action?pageId=78644755&src=contextnavpagetreemode)

*   [Support and Contacts](https://wiki.nci.nih.gov/display/TCGA/Support+and+Contacts?src=contextnavpagetreemode)






[](https://wiki.nci.nih.gov/collector/pages.action?key=TCGA&src=sidebar " (g , s)")<h2>Space Details</h2><div class="personal-space-logo-hint">Your profile picture is used as the logo for your personal space. <a href="/users/profile/editmyprofilepicture.action" target="_blank">Change your profile picture</a>.</div>

Reorder pages

ConfigureSpace tools

*   [Overview](https://wiki.nci.nih.gov/spaces/viewspacesummary.action?key=TCGA&src=spacetools)
*   [Content Tools](https://wiki.nci.nih.gov/pages/reorderpages.action?key=TCGA&src=spacetools)

*   [Reorder pages](https://wiki.nci.nih.gov/pages/reorderpages.action?key=TCGA&src=spacetools)

*   [Save for later](https://wiki.nci.nih.gov/labels/addfavourite.action?entityId=26187384&atl_token=3727b31847763e4f3b4f6c758b5542dc19c12d84 "Save for later (f)")
*   [Watch](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "Watch (w)")
*   [Share](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "Share this content by emailing other users (s)")
*   [](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   *   [Attachments (2)](https://wiki.nci.nih.gov/pages/viewpageattachments.action?pageId=26187384 "View Attachments (t)")
   *   [Page History](https://wiki.nci.nih.gov/pages/viewpreviousversions.action?pageId=26187384)
   *   [Scaffolding History Beta](https://wiki.nci.nih.gov/pages/viewscaffoldversions.action?pageId=26187384)
   *   [Restrictions](https://wiki.nci.nih.gov/pages/viewinfo.action?pageId=26187384 "Edit restrictions")

   *   [Page Information](https://wiki.nci.nih.gov/pages/viewinfo.action?pageId=26187384)
   *   [Resolved comments (0)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Link to this Page…](https://wiki.nci.nih.gov/pages/viewinfo.action?pageId=26187384 "Link to this Page (k)")
   *   [View in Hierarchy](https://wiki.nci.nih.gov/pages/reorderpages.action?key=TCGA&openId=26187384#selectedPageInHierarchy)
   *   [View Storage Format](https://wiki.nci.nih.gov/plugins/viewstorage/viewpagestorage.action?pageId=26187384)
   *   [View Source](https://wiki.nci.nih.gov/plugins/viewsource/viewpagesrc.action?pageId=26187384)
   *   [View Scaffolding XML](https://wiki.nci.nih.gov/pages/metadata/viewxml.action?pageId=26187384)
   *   [Export to PDF](https://wiki.nci.nih.gov/spaces/flyingpdf/pdfpageexport.action?pageId=26187384)
   *   [Export to Word](https://wiki.nci.nih.gov/exportword?pageId=26187384)
   *   [View Visio File](https://wiki.nci.nih.gov/plugins/lucidchart/selectVisio.action?contentId=26187384)

   *   [Copy Page Tree](https://wiki.nci.nih.gov/plugins/tree-copy/preparing-copying.action?pageId=26187384)
   *   [Delete Page Tree](https://wiki.nci.nih.gov/plugins/tree-delete/preparing-delete.action?pageId=26187384)
   *   [Copy](https://wiki.nci.nih.gov/pages/copypage.action?idOfPageToCopy=26187384&spaceKey=TCGA)
   *   [Copy with Scaffolding XML](https://wiki.nci.nih.gov/pages/copyscaffoldfromajax.action?idOfPageToCopy=26187384&spaceKey=TCGA)


1.  [Pages](https://wiki.nci.nih.gov/collector/pages.action?key=TCGA&src=breadcrumbs-collector)
2.  **…**
3.  [TCGA Wiki Home](https://wiki.nci.nih.gov/display/TCGA/TCGA+Wiki+Home?src=breadcrumbs-expanded)
4.  [TCGA User Documentation](https://wiki.nci.nih.gov/display/TCGA/TCGA+User+Documentation?src=breadcrumbs-expanded)
5.  [TCGA Specifications](https://wiki.nci.nih.gov/display/TCGA/TCGA+Specifications?src=breadcrumbs-expanded)
6.  [File Format Specifications](https://wiki.nci.nih.gov/display/TCGA/File+Format+Specifications?src=breadcrumbs-expanded)
7.  [Mutation Annotation Format (MAF) Specification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification?src=breadcrumbs-parent)

[Skip to end of banner](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#page-banner-end)

*   [JIRA links](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4)

[Go to start of banner](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#page-banner-start)

[Mutation Annotation Format (MAF) Specification - v2.4](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4)
========================================================================================================================================================

<table class="aui"> <thead> <tr class="header"> <th class="search-result-title">Page Title</th> <th class="search-result-space">Space</th> <th class="search-result-date">Updated</th> </tr> </thead> </table> <p class="search-result-count">{0}</p> <tr class="search-result"> <td class="search-result-title"><a href="{1}" class="content-type-{2}"><span>{0}</span></a></td> <td class="search-result-space"><a class="space" href="/display/{4}/" title="{3}">{3}</a></td> <td class="search-result-date"><span class="date" title="{6}">{5}</span></td> </tr> [Skip to end of metadata](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#page-metadata-end)

*   Created by Unknown User (arik), last modified by Unknown User (ayalabe) on [Jul 09, 2014](https://wiki.nci.nih.gov/pages/diffpagesbyversion.action?pageId=26187384&selectedPageVersions=96&selectedPageVersions=97 "Show changes")

[Go to start of metadata](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#page-metadata-start)

**IMPORTANT**: MAF files can be submitted to the DCC **only** by following the procedure described [here](https://wiki.nci.nih.gov/x/yw7RB).

**Document Information**

**Specification for Mutation Annotation Format**  
Version 2.4  
March 6, 2013

**Contents**

*   1[Current version changes](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Currentversionchanges)
*   2[About MAF specifications](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-AboutMAFspecifications)
   *   2.1[Definition of open access MAF data](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-DefinitionofopenaccessMAFdata)
   *   2.2[Somatic MAF vs. Protected MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-SomaticMAFvs.ProtectedMAF)
*   3[MAF file fields](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-MAFfilefields)
   *   3.1[Table 1 - File column headers](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Table1-Filecolumnheaders)
   *   3.2[MAF file checks](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-MAFfilechecks)
*   4[MAF naming convention](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-MAFnamingconvention)
*   5[Previous specification versions](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Previousspecificationversions)

Current version changes
=======================

This current revision is **version 2.4** of the Mutation Annotation Format (MAF) specification.

The following items in the specification were added or modified in version 2.4 from version 2.3:

*   Header for MAF file is "#version 2.4"
*   Allowed "[Validation_Status](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Validation_Status)" values are now "Untested, Inconclusive, Valid, and Invalid."
*   "[Mutation_Status](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Mutation_Status)" values are now aligned with the [VCF VLS](https://wiki.nci.nih.gov/x/2gcYAw) field and are "None, Germline, Somatic, LOH, Post-transcriptional modification, Unknown" (addition of the enumerated value of "Post-transcriptional modification")  
   *   The VCF VLS value of Wild Type will be replaced by None
*   The value in the "[Validation_Status](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Validation_Status)" field determines what values are allowed in the "[Mutation_Status](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Mutation_Status)" field
*   The "[Validation_Status](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Validation_Status)" field will no longer accept NULL values
*   The values allowed in the "[Sequence_Source](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Sequence_Source)" column have been changed to be a subset of the [SRA 1.5 library_strategy values](http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.experiment.xsd?view=markup).
*   All Version 2.4 MAF files will be required to contain UUID columns (Tumor\_Sample\_UUID and Matched\_Norm\_Sample_UUID)
*   For the "[Variant_Classification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-Variant_Classification)" field, the values of **De\_novo\_Start_InFrame** and **De\_novo\_Start_OutOfFrame** are no longer allowed.
*   The values for the"[dbSNP\_Val\_Status](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-dbSNP_Val_Status)" will be enforced. "none" will no longer be allowed.
*   "Somatic" is the only acceptable value for 'Mutation\_Status' for a somatic MAF (named .somatic.maf). Protected MAF (named .protected.maf) has no such restriction and can contain Somatic, Germline, Unknown, LOH, Post-transcriptional modification, None for Mutation\_Status.
*   For a somatic MAF, following rules should be satisfied:  
   SOMATIC = A AND (B OR C OR D)  
   A: _Mutation_Status_ == "Somatic"  
   B: _Validation_Status_ == "Valid"  
   C. _Verification_Status_ == "Verified"  
   D. _Variant_Classification_ is not {Intron, 5'UTR, 3'UTR, 5'Flank, 3'Flank, IGR}, which implies that _Variant_Classification_ can only be \\{Frame\_Shift\_Del, Frame\_Shift\_Ins, In\_Frame\_Del, In\_Frame\_Ins, Missense\_Mutation, Nonsense\_Mutation, Silent, Splice\_Site, Translation\_Start\_Site, Nonstop\_Mutation, RNA, Targeted_Region}.
*   If Validation\_Status  == "Valid" then Tumor\_Validation\_Allele1, Tumor\_Validation\_Allele2, Match\_Norm\_Validation\_Allele1, Match\_Norm\_Validation_Allele2 cannot  be null 
*   MAF 2.3 Validation rule 8a was superceded by:  
   *   **If Validation\_Status == "Invalid" then Tumor\_Validation\_Allele1, Tumor\_Validation\_Allele2, Match\_Norm\_Validation\_Allele1, Match\_Norm\_Validation_Allele2 cannot be null AND Tumor\_Validation\_Allelle1 == Match\_Norm\_Validation\_Allele1 AND Tumor\_Validation\_Allelle2 == Match\_Norm\_Validation\_Allele2  (Added as a replacement for 8a as a result of breakdown)**  
          



About MAF specifications
========================

Mutation annotation files should be transferred to the DCC. Those files should be formatted using the mutation annotation format (MAF) that is described below. File naming convention is also [below](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#MutationAnnotationFormat(MAF)Specification-v2.4-MAFnamingconvention).

Following categories of somatic mutations are reported in MAF files:

*   Missense and nonsense
*   Splice site, defined as SNP within 2 bp of the splice junction
*   Silent mutations
*   Indels that overlap the coding region or splice site of a gene or the targeted region of a genetic element of interest.
*   Frameshift mutations
*   Mutations in regulatory regions

### Definition of open access MAF data

A large proportion of MAFs are submitted as discovery data and sites labeled as somatic in these files overlap with known germline variants. In order to minimize germline contamination in putative (unvalidated) somatic calls, certain filtering criteria have been imposed. Based on current policy, open access MAF data should:

*   **include** all validated somatic mutation calls
*   **include** all unvalidated somatic mutation calls that overlap with a coding region or splice site
*   **exclude** all other types of mutation calls (i.e., non-somatic calls (validated or not), unvalidated somatic calls that are not in coding region or splice sites, and dbSNP sites that are not annotated as somatic in dbSNP, COSMIC or OMIM)

### Somatic MAF vs. Protected MAF

Centers will submit to the DCC MAF archives that contain Somatic MAF (named **.somatic.maf**) for open access data and an all-inclusive Protected MAF (named **.protected.maf**) that does not filter any data out and represents the original super-set of mutation calls. The files will be formatted using the Mutation Annotation Format (MAF).

The following table lists some of the critical attributes of somatic and protected MAF files and provides a comparison.

**Attribute**

Somatic MAF

Protected MAF

**Attribute**

Somatic MAF

Protected MAF

**File naming**

Somatic MAFs should be named as ***.somatic.maf** and cannot contain 'germ' or 'protected' in file name.

Protected MAFs should be named as***.protected.maf** and should not contain 'somatic' in the file name.

**Mutation category**

Somatic MAFs can only contain entries where _Mutation_Status_ is "Somatic". If any other value is assigned to the field, the archive will fail. Experimentally validated or unvalidated (see next row) somatic mutations can be included in the file.

There is no such restriction for protected MAF. The file should contain all mutation calls including those from which .somatic.maf is derived.

**Filtering criteria**

In order to minimize germline contamination, somatic MAFs can contain unvalidated somatic mutations only from coding regions and splice sites, which implies:   
If _Validation_Status_ **is** "Unknown", _V__a__riant_Classification_ **cannot** be 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR, or Intron. _Variant_Classification_ can only be \\{Frame\_Shift\_Del, Frame\_Shift\_Ins, In\_Frame\_Del, In\_Frame\_Ins, Missense\_Mutation, Nonsense\_Mutation, Silent, Splice\_Site, Translation\_Start\_Site, Nonstop\_Mutation, RNA, Targeted\_Region, De\_novo\_Start\_InFrame, De\_novo\_Start_OutOfFrame\\}.  
There is no such constraint for experimentally validated (_Validation_Status_ is "Valid") somatic mutations.   

dbSNP sites that are not annotated as somatic in dbSNP, COSMIC or OMIM must be removed from somatic MAFs.

There are no such constraints for mutations in protected MAF.

**Access level**

These files are deployed as open access data.

These files are deployed as protected data.

MAF file fields
===============

The format of a MAF file is tab-delimited columns. Those columns are described in Table 1 and are required in every MAF file. The order of the columns will be validated by the DCC. Column headers and values **are** case sensitive where specified. Columns may allow null values (i.e._ blank cells) and/or have enumerated values. **The validator looks for a header stating the version of the specification to validate against (e.g. #version 2.4). If not, validation fails.** Any columns that come after the columns described in Table 1 are optional. Optional columns are not validated by the DCC and can be in any order.

Table 1 - File column headers
-----------------------------

**Index**

**MAF Column Header**

**Description of Values**


**Example**

**Case**  
**Sensitive**

**Null**

**Enumerated**

1

Hugo_Symbol

HUGO symbol for the gene (HUGO symbols are _always_ in all caps). If no gene exists within 3kb enter "Unknown".  
Source: [http://genenames.org](http://genenames.org/)

EGFR

Yes

No

Set or Unknown

2

Entrez\_Gene\_Id

Entrez gene ID (an integer). If no gene exists within 3kb enter "0".  
Source: [http://ncbi.nlm.nih.gov/sites/entrez?db=gene](http://ncbi.nlm.nih.gov/sites/entrez?db=gene)

1956

No

No

Set

3

Center

Genome sequencing center reporting the variant. If multiple institutions report the same mutation separate list using semicolons. Non-GSC centers will be also supported if center name is an accepted center name.

hgsc.bcm.edu;genome.wustl.edu

Yes

No

Set

4

NCBI_Build

Any TGCA accepted genome identifier.  Can be string, integer or a float.  

hg18, hg19, GRCh37, GRCh37-lite, 36, 36.1, 37,


No

No

Set and Enumerated.

5

Chromosome

Chromosome number without "chr" prefix that contains the gene.

X, Y, M, 1, 2, etc.

Yes

No

Set

6

Start_Position

Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate (1-based coordinate system).

999

No

No

Set

7

End_Position

Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate (inclusive, 1-based coordinate system).

1000

No

No

Set

8

Strand

Genomic strand of the reported allele. Variants should always be reported on the positive genomic strand. (Currently, only the positive strand is an accepted value).

+

No

No

+

9

Variant_Classification

Translational effect of variant allele.

Missense_Mutation

Yes

No

Frame\_Shift\_Del, Frame\_Shift\_Ins, In\_Frame\_Del, In\_Frame\_Ins, Missense\_Mutation, Nonsense\_Mutation, Silent, Splice\_Site, Translation\_Start\_Site, Nonstop\_Mutation, 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR 1 , Intron, RNA, Targeted_Region

10

Variant_Type

Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP but for 3 consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of 4 or more.

INS

Yes

No

SNP, DNP, TNP, ONP, INS, DEL, or Consolidated 2

11

Reference_Allele

The plus strand reference allele at this position. Include the sequence deleted for a deletion, or "-" for an insertion.

A

Yes

No

A,C,G,T and/or -

12

Tumor\_Seq\_Allele1

Primary data genotype. Tumor sequencing (discovery) allele 1. " -" for a deletion represent a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.

C

Yes

No

A,C,G,T and/or -

13

Tumor\_Seq\_Allele2

Primary data genotype. Tumor sequencing (discovery) allele 2. " -" for a deletion represents a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.

G

Yes

No

A,C,G,T and/or -

14

dbSNP_RS

Latest dbSNP rs ID (dbSNP_ID) or "novel" if there is no dbSNP record. source: [http://ncbi.nlm.nih.gov/projects/SNP/](http://ncbi.nlm.nih.gov/projects/SNP/)

rs12345

Yes

Yes

Set or "novel"

15

dbSNP\_Val\_Status

dbSNP validation status. Semicolon- separated list of validation statuses.

by2Hit2Allele;byCluster

No

Yes

by1000genomes;by2Hit2Allele; byCluster; byFrequency; byHapMap; byOtherPop; bySubmitter; alternate_allele 3 **Note that "none" will no longer be an acceptable value.**

16

Tumor\_Sample\_Barcode

BCR aliquot barcode for the tumor sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID.

TCGA-02-0021-01A-01D-0002-04

Yes

No

Set

17

Matched\_Norm\_Sample_Barcode

BCR aliquot barcode for the matched normal sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID; e.g. TCGA-02-0021-10A-01D-0002-04 (compare portion ID '10A' normal sample, to '01A' tumor sample).

TCGA-02-0021-10A-01D-0002-04

Yes

No

Set

18

Match\_Norm\_Seq_Allele1

Primary data. Matched normal sequencing allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.

T

Yes

Yes

A,C,G,T and/or -

19

Match\_Norm\_Seq_Allele2

Primary data. Matched normal sequencing allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.

ACGT

Yes

Yes

A,C,G,T and/or -

20

Tumor\_Validation\_Allele1

Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.

-

Yes

Yes

A,C,G,T and/or -

21

Tumor\_Validation\_Allele2

Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.

A

Yes

Yes

A,C,G,T and/or -

22

Match\_Norm\_Validation_Allele1

Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.

C

Yes

Yes

A,C,G,T and/or -

23

Match\_Norm\_Validation_Allele2

Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.

G

Yes

Yes

A,C,G,T and/or -

24

Verification_Status 4

Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing.



Verified

Yes

Yes

Verified, Unknown

25

Validation_Status  5

Second pass results from orthogonal technology.

Valid

Yes

No

*   Untested
*   Inconclusive
*   Valid
*   Invalid

26

Mutation_Status

Updated to reflect validation or verification status and to be in agreement with the [VCF VLS](https://wiki.nci.nih.gov/x/2gcYAw) field. The values allowed in this field are constrained by the value in the Validation_Status field.

Somatic

Yes

No

Validation_Status value

Allowed 6 Mutation_Status values

Untested

*   None
*   Germline
*   Somatic
*   LOH
*   Post-transcriptional modification
*   Unknown

Inconclusive

*   None
*   Germline
*   Somatic
*   LOH
*   Post-transcriptional modification
*   Unknown

Valid

*   Germline
*   Somatic
*   LOH
*   Post-transcriptional modification
*   Unknown

Invalid

*   None  



27

Sequencing_Phase

TCGA sequencing phase. Phase should change under any circumstance that the targets under consideration change.

Phase_I

No

Yes

No

28

Sequence_Source

Molecular assay type used to produce the analytes used for sequencing. Allowed values are a subset of the [SRA 1.5](http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/) library_strategy field values. This subset matches those used at CGHub.

WGS;WXS

Yes

No

*   Common TCGA values  
   *   WGS
   *   WGA
   *   WXS
   *   RNA-Seq
   *   miRNA-Seq
   *   Bisulfite-Seq
   *   VALIDATION
   *   Other
*   Other allowed values (per SRA 1.5)  
   *   ncRNA-Seq
   *   WCS
   *   CLONE
   *   POOLCLONE
   *   AMPLICON
   *   CLONEEND
   *   FINISHING
   *   ChIP-Seq
   *   MNase-Seq
   *   DNase-Hypersensitivity
   *   EST
   *   FL-cDNA
   *   CTS
   *   MRE-Seq
   *   MeDIP-Seq
   *   MBD-Seq
   *   Tn-Seq
   *   FAIRE-seq
   *   SELEX
   *   RIP-Seq
   *   ChIA-PET

29

Validation_Method

The assay platforms used for the validation call. Examples: Sanger\_PCR\_WGA, Sanger\_PCR\_gDNA, 454\_PCR\_WGA, 454\_PCR\_gDNA; separate multiple entries using semicolons.

Sanger\_PCR\_WGA;Sanger\_PCR\_gDNA

No

**NO**. I**f Validation_Status = Untested then "none"**

No

30

Score

Not in use.

NA

No

Yes

No

31

BAM_File

Not in use.

NA

No

Yes

No

32

Sequencer

Instrument used to produce primary data. Separate multiple entries using semicolons.

Illumina GAIIx;SOLID

Yes

No

*   Illumina GAIIx
*   Illumina HiSeq
*   SOLID
*   454
*   ABI 3730xl
*   Ion Torrent PGM
*   Ion Torrent Proton
*   PacBio RS
*   Illumina MiSeq
*   Illumina HiSeq 2500
*   454 GS FLX Titanium
*   AB SOLiD 4 System

33

Tumor\_Sample\_UUID

BCR aliquot UUID for tumor sample

550e8400-e29b-41d4-a716-446655440000

Yes

No

 

34

Matched\_Norm\_Sample_UUID

BCR aliquot UUID for matched normal

567e8487-e29b-32d4-a716-446655443246

Yes

No

 

1

Intergenic Region

↩  
2

'Consolidated' is used to indicate a site that was initially reported as as variant but subsequently removed from further analysis because it was consolidated into a new variant.  For example, a SNP variant incorporated into a TNP variant.

↩  
3

Used when the discovered variant differs from that of dbSNP

↩  
4

These MAF headers describe the technology that was used to confirm a mutation, whether the same technology (“verification”) or a different technology (“validation”) is used to prove that a variant is germline or a somatic mutation.

↩  
5

These MAF headers describe the technology that was used to confirm a mutation, whether the same technology (“verification”) or a different technology (“validation”) is used to prove that a variant is germline or a somatic mutation.

↩  
6

Explanation of some Validation Status-Mutation Status combinations

Validation Status

Mutation Status

Explanation

Validation Status

Mutation Status

Explanation

Valid

Unknown

a valid variant with unknown somatic status due to lack of data from matched normal tissue.

Invalid

None

validation attempted, tumor and normal are homozygous reference (formerly described as Wildtype)

Inconclusive

Unknown

validation failed, neither the genotype nor its somatic status is certain due to lack of data from matched normal tissue

Inconclusive

None

validation failed, tumor genotype appears to be homozygous reference

↩  

Important Criteria

**Index column indicates the order in which the columns are expected**. **All headers are case sensitive.** The Case Sensitive column specifies which values are case sensitive. The Null column indicates which MAF columns are allowed to have null values. The Enumerated column indicates which MAF columns have specified values: an Enumerated value of "No" indicates that there are no specified values for that column; other values indicate the specific values listed allowed; a value of "Set" indicates that the MAF column values come from a specified set of known values (_e.g._ HUGO gene symbols).

MAF file checks
---------------

The DCC Archive Validator checks the integrity of a MAF file. Validation will fail if any of the below are not true for a MAF file:

1.  Column header text (including case) and order must match specification (Table 1) exactly
2.  Values under column headers listed in the specification (Table 1) as not null must have values
3.  Values that are specified in Table 1 as Case Sensitive must be.
4.  If column headers are listed in the specification as having _enumerated_ values (_i.e._ a "Yes" in the "Enumerated" column), then the values under those column must come from the enumerated values listed under "Enumerated".
5.  If column headers are listed in the specification as having _set_ values (_i.e._ a "Set" in the "Enumerated" column), then the values under those column must come from the enumerated values of that domain (_e.g._ HUGO gene symbols).
6.  All Allele-based columns must contain - (deletion), or a string composed of the following capitalized letters: A, T, G, C.
7.  If Validation\_Status == "Untested" then Tumor\_Validation\_Allele1, Tumor\_Validation\_Allele2, Match\_Norm\_Validation\_Allele1, Match\_Norm\_Validation\_Allele2 can be null (depending on Validation\_Status).  

   1.  If Validation\_Status == "Inconclusive" then Tumor\_Validation\_Allele1, Tumor\_Validation\_Allele2, Match\_Norm\_Validation\_Allele1, Match\_Norm\_Validation\_Allele2 can be null (depending on Validation\_Status)**.**
8.  If Validation\_Status == Valid, then Validated\_Tumor\_Allele1 and Validated\_Tumor_Allele2must be populated (one of A, C, G, T, and -)  
   1.  If Validation\_Status  == "Valid" then Tumor\_Validation\_Allele1, Tumor\_Validation\_Allele2, Match\_Norm\_Validation\_Allele1, Match\_Norm\_Validation_Allele2 cannot  be null 
   2.   If Validation\_Status == "Invalid" then Tumor\_Validation\_Allele1, Tumor\_Validation\_Allele2, Match\_Norm\_Validation\_Allele1, Match\_Norm\_Validation_Allele2 cannot be null AND Tumor\_Validation\_Allelle1 == Match\_Norm\_Validation\_Allele1 AND Tumor\_Validation\_Allelle2 == Match\_Norm\_Validation\_Allele2  (Added as a replacement for 8a as a result of breakdown)
9.  Check allele values against Mutation_Status:  
   Check allele values against Validation_status:
   1.  If Mutation_Status == "Germline" and Validation_Status == "Valid", then Tumor\_Validation\_Allele1 == Match\_Norm\_Validation\_Allele1 and Tumor\_Validation\_Allele2 == Match\_Norm\_Validation\_Allele2.  


   2.  If Mutation\_Status == "Somatic" and Validation\_Status == "Valid", then Match\_Norm\_Validation\_Allele1 == Match\_Norm\_Validation\_Allele2 == Reference\_Allele and (Tumor\_Validation\_Allele1 or Tumor\_Validation\_Allele2) != Reference\_Allele  


   3.  If Mutation\_Status == "LOH" and Validation\_Status=="Valid", then Tumor\_Validation\_Allele1 == Tumor\_Validation\_Allele2 and Match\_Norm\_Validation\_Allele1 != Match\_Norm\_Validation\_Allele2 and Tumor\_Validation\_Allele1 == (Match\_Norm\_Validation\_Allele1 or Match\_Norm\_Validation\_Allele2).  


10.  Check that Start\_position <= End\_position
11.  Check for the Start\_position and End\_position against Variant_Type:  
   1.  If Variant\_Type == "INS", then (End\_position - Start\_position + 1 == length (Reference\_Allele) or End\_position - Start\_position == 1) and length(Reference\_Allele) <= length(Tumor\_Seq\_Allele1 and Tumor\_Seq_Allele2)
   2.  If Variant\_Type == "DEL", then End\_position - Start\_position + 1 == length (Reference\_Allele), then length(Reference\_Allele) >= length(Tumor\_Seq\_Allele1 and Tumor\_Seq_Allele2)
   3.  If Variant\_Type == "SNP", then length(Reference\_Allele and Tumor\_Seq\_Allele1 and Tumor\_Seq\_Allele2) ==  1 and (Reference\_Allele and Tumor\_Seq\_Allele1 and Tumor\_Seq_Allele2) != "-"
   4.  If Variant\_Type == "DNP", then length(Reference\_Allele and Tumor\_Seq\_Allele1 and Tumor\_Seq\_Allele2) ==  2 and (Reference\_Allele and Tumor\_Seq\_Allele1 and Tumor\_Seq_Allele2) !contain "-"
   5.  If Variant\_Type == "TNP", then length(Reference\_Allele and Tumor\_Seq\_Allele1 and Tumor\_Seq\_Allele2) ==  3 and (Reference\_Allele and Tumor\_Seq\_Allele1 and Tumor\_Seq_Allele2) !contain "-"
   6.  If Variant\_Type == "ONP", then length(Reference\_Allele) == length(Tumor\_Seq\_Allele1) == length(Tumor\_Seq\_Allele2) > 3 and (Reference\_Allele and Tumor\_Seq\_Allele1 and Tumor\_Seq_Allele2) !contain "-"
12.  Validation for UUID-based files:  
   1.  Column #33 must be Tumor\_Sample\_UUID containing UUID of the BCR aliquot for tumor sample
   2.  Column #34 must be Matched\_Norm\_Sample_UUID containing UUID of the BCR aliquot for matched normal sample
   3.  Metadata represented by Tumor\_Sample\_Barcode and Matched\_Norm\_Sample\_Barcode should correspond to the UUIDs assigned to Tumor\_Sample\_UUID and Matched\_Norm\_Sample\_UUID respectively

MAF naming convention
=====================

In archives uploaded to the DCC, the MAF file name should relate to the containing archive name in the following way:

If the archive has the name

[?](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

`<domain>_<disease_abbrev>.<platform>.Level_2.<serial_index>.<revision>.0.tar.gz`

then a somatic MAF file with the archive should be named according to

[?](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

`<domain>_<disease_abbrev>.<platform>.Level_2.<serial_index>[.<optional_tag>].somatic.maf`

and a protected MAF with the archive should be named according to

[?](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

`<domain>_<disease_abbrev>.<platform>.Level_2.<serial_index>[.<optional_tag>].``protected``.maf`

The `<optional_tag>` may consist of alphanumeric characters, dash, and underscore; no spaces or periods; or it may be left out altogether. The purpose of the optional tag is to impart some brief annotation.

_Example_

For the archive

[?](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

`genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.7.6.0.tar.gz`

the following are examples of valid maf names

[?](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

`genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.7.preliminary.somatic.maf`

`genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.7.protected.maf`

Previous specification versions
===============================

*   [version 2.3](https://wiki.nci.nih.gov/x/xYGDBw)
*   [version 2.1](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.2)
*   [version 2.1](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.1)
*   [version 2.0](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.0)
*   [version 1.0](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v1.0)

*   [binf](https://wiki.nci.nih.gov/label/TCGA/binf)
*   [latest](https://wiki.nci.nih.gov/label/TCGA/latest)
*   [Edit Labels](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4# "Edit Labels (l)")

[![User icon: Add a picture of yourself](./Mutation Annotation Format (MAF) Specification - v2.4 - TCGA - National Cancer Institute - Confluence Wiki_files/add_profile_pic.png)](https://wiki.nci.nih.gov/users/profile/editmyprofilepicture.action)

Write a comment…

[Add Comment](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4?showComments=true&showCommentArea=true#addcomment)

Overview

Content Tools

*   Powered by [Atlassian Confluence](http://www.atlassian.com/software/confluence) 5.10.4
*   Printed by Atlassian Confluence 5.10.4

[Atlassian](http://www.atlassian.com/)

var \_gaq = \_gaq || \[\]; \_gaq.push(\['\_setAccount', 'UA-4968072-3'\]); \_gaq.push(\['\_trackPageview'\]); (function() { var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true; ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js'; var s = document.getElementsByTagName('script')\[0\]; s.parentNode.insertBefore(ga, s); })();

{"serverDuration": 531, "requestCorrelationId": "b63c2888df75e7ca"} AJS.BigPipe = AJS.BigPipe || {}; AJS.BigPipe.metrics = AJS.BigPipe.metrics || {}; AJS.BigPipe.metrics.pageEnd = typeof window.performance !== "undefined" && typeof window.performance.now === "function" ? Math.ceil(window.performance.now()) : 0; AJS.BigPipe.metrics.isBigPipeEnabled = 'false' === 'true';

1.  [TCGA](https://wiki.nci.nih.gov/display/TCGA)
2.  [Pages](https://wiki.nci.nih.gov/collector/pages.action?key=TCGA)
3.  **…**
4.  [TCGA Wiki Home](https://wiki.nci.nih.gov/display/TCGA/TCGA+Wiki+Home)
5.  [TCGA User Documentation](https://wiki.nci.nih.gov/display/TCGA/TCGA+User+Documentation)
6.  [TCGA Specifications](https://wiki.nci.nih.gov/display/TCGA/TCGA+Specifications)
7.  [File Format Specifications](https://wiki.nci.nih.gov/display/TCGA/File+Format+Specifications)
8.  [Mutation Annotation Format (MAF) Specification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification)
9.  [Mutation Annotation Format (MAF) Specification - v2.4](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4)

<link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/1.1.11/\_/download/batch/com.atlassian.confluence.plugins.confluence-collaborative-editor-plugin:confluence-collaborative-editor-plugin-editor-content-resources/com.atlassian.confluence.plugins.confluence-collaborative-editor-plugin:confluence-collaborative-editor-plugin-editor-content-resources.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-collaborative-editor-plugin:confluence-collaborative-editor-plugin-editor-content-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/4.1.3/\_/download/batch/com.stiltsoft.confluence.plugin.tablefilter.tablefilter:tf-editor-content-resources/com.stiltsoft.confluence.plugin.tablefilter.tablefilter:tf-editor-content-resources.css" data-wrm-key="com.stiltsoft.confluence.plugin.tablefilter.tablefilter:tf-editor-content-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/a191f6d5283eab057f5eb0b18f7e7453-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.9.24/\_/download/batch/com.atlassian.auiplugin:aui-page-typography/com.atlassian.auiplugin:aui-page-typography.css" data-wrm-key="com.atlassian.auiplugin:aui-page-typography" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/a191f6d5283eab057f5eb0b18f7e7453-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.9.24/\_/download/batch/com.atlassian.auiplugin:aui-avatars/com.atlassian.auiplugin:aui-avatars.css" data-wrm-key="com.atlassian.auiplugin:aui-avatars" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/202a2990a81c67383e0548a13060e251-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.9.24/\_/download/batch/com.atlassian.auiplugin:aui-page-layout/com.atlassian.auiplugin:aui-page-layout.css" data-wrm-key="com.atlassian.auiplugin:aui-page-layout" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/0ff1e2a127cdaebdc2dd0bdef1cf04fe-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.10.4/\_/download/batch/com.atlassian.confluence.editor:editor-content-styles/com.atlassian.confluence.editor:editor-content-styles.css" data-wrm-key="com.atlassian.confluence.editor:editor-content-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/0ff1e2a127cdaebdc2dd0bdef1cf04fe-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.10.4/\_/download/batch/com.atlassian.confluence.editor:table-resizable-editor-content-styles/com.atlassian.confluence.editor:table-resizable-editor-content-styles.css?confluence.table.resizable=true" data-wrm-key="com.atlassian.confluence.editor:table-resizable-editor-content-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/0ff1e2a127cdaebdc2dd0bdef1cf04fe-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.10.4/\_/download/batch/com.atlassian.confluence.editor:table-resizable-styles/com.atlassian.confluence.editor:table-resizable-styles.css?confluence.table.resizable=true" data-wrm-key="com.atlassian.confluence.editor:table-resizable-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.10.4/\_/download/batch/com.atlassian.confluence.plugins.confluence-templates:variable-editor-content-styles/com.atlassian.confluence.plugins.confluence-templates:variable-editor-content-styles.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-templates:variable-editor-content-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/a191f6d5283eab057f5eb0b18f7e7453-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.9.24/\_/download/batch/com.atlassian.auiplugin:aui-lozenge/com.atlassian.auiplugin:aui-lozenge.css" data-wrm-key="com.atlassian.auiplugin:aui-lozenge" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/2.17/\_/download/batch/com.atlassian.confluence.plugins.status-macro:view\_content\_status/com.atlassian.confluence.plugins.status-macro:view\_content\_status.css" data-wrm-key="com.atlassian.confluence.plugins.status-macro:view\_content\_status" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/489901e715cade6fe0efc7e21e194576-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/2.17/\_/download/batch/com.atlassian.confluence.plugins.status-macro:editor\_content\_status/com.atlassian.confluence.plugins.status-macro:editor\_content\_status.css" data-wrm-key="com.atlassian.confluence.plugins.status-macro:editor\_content\_status" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d0638b94054c44b2ca8199723136b528-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/1.2.8/\_/download/batch/com.atlassian.confluence.plugins.confluence-fixed-headers:confluence-fixed-headers-editor-content-resources/com.atlassian.confluence.plugins.confluence-fixed-headers:confluence-fixed-headers-editor-content-resources.css?confluence.view.edit.transition=true" data-wrm-key="com.atlassian.confluence.plugins.confluence-fixed-headers:confluence-fixed-headers-editor-content-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/5baeb58a5608e28c5f241eee74f064b2-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/4.0.2/\_/download/batch/confluence.extra.attachments:attachments-css/confluence.extra.attachments:attachments-css.css" data-wrm-key="confluence.extra.attachments:attachments-css" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/2.6.1/\_/download/batch/nl.avisi.confluence.plugins.numberedheadings:nh-tinymce-css-resources/nl.avisi.confluence.plugins.numberedheadings:nh-tinymce-css-resources.css" data-wrm-key="nl.avisi.confluence.plugins.numberedheadings:nh-tinymce-css-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/13.0.1/\_/download/batch/com.atlassian.confluence.plugins.confluence-roadmap-plugin:roadmap-placeholder-css/com.atlassian.confluence.plugins.confluence-roadmap-plugin:roadmap-placeholder-css.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-roadmap-plugin:roadmap-placeholder-css" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/1.0/\_/download/batch/confluence.web.resources:panel-styles/confluence.web.resources:panel-styles.css" data-wrm-key="confluence.web.resources:panel-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/66f9fcaa7e13c424bc045a814fa5aefe-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/1.0/\_/download/batch/confluence.web.resources:content-styles/confluence.web.resources:content-styles.css" data-wrm-key="confluence.web.resources:content-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/1.23.10/\_/download/batch/com.lucidchart.confluence.plugins.lucid-confluence:lucidchart-editor-resources/com.lucidchart.confluence.plugins.lucid-confluence:lucidchart-editor-resources.css" data-wrm-key="com.lucidchart.confluence.plugins.lucid-confluence:lucidchart-editor-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/8.17.8/\_/download/batch/net.customware.confluence.plugin.scaffolding:editor-content-resources/net.customware.confluence.plugin.scaffolding:editor-content-resources.css" data-wrm-key="net.customware.confluence.plugin.scaffolding:editor-content-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.10.4/\_/download/batch/com.atlassian.confluence.plugins.confluence-page-layout:editor-pagelayout-content-styles/com.atlassian.confluence.plugins.confluence-page-layout:editor-pagelayout-content-styles.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-page-layout:editor-pagelayout-content-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/4.1.4/\_/download/batch/com.atlassian.confluence.plugins.confluence-business-blueprints:sharelinks-urlmacro-editor-resources/com.atlassian.confluence.plugins.confluence-business-blueprints:sharelinks-urlmacro-editor-resources.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-business-blueprints:sharelinks-urlmacro-editor-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/7.0.2/\_/download/batch/com.atlassian.confluence.plugins.confluence-inline-tasks:inline-tasks-styles/com.atlassian.confluence.plugins.confluence-inline-tasks:inline-tasks-styles.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-inline-tasks:inline-tasks-styles" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/7.0.2/\_/download/batch/com.atlassian.confluence.plugins.confluence-inline-tasks:editor-autocomplete-date-css/com.atlassian.confluence.plugins.confluence-inline-tasks:editor-autocomplete-date-css.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-inline-tasks:editor-autocomplete-date-css" data-wrm-batch-type="resource" media="all"> <!--\[if lte IE 9\]> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/7.0.2/\_/download/batch/com.atlassian.confluence.plugins.confluence-inline-tasks:editor-autocomplete-date-css/com.atlassian.confluence.plugins.confluence-inline-tasks:editor-autocomplete-date-css.css?conditionalComment=lte+IE+9" data-wrm-key="com.atlassian.confluence.plugins.confluence-inline-tasks:editor-autocomplete-date-css" data-wrm-batch-type="resource" media="all"> <!\[endif\]--> <link type="text/css" rel="stylesheet" href="/s/51a3950c50505d7f5d8217ed8913f870-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/1.5.5/\_/download/batch/com.atlassian.confluence.plugins.confluence-view-file-macro:view-file-macro-editor-content-resources/com.atlassian.confluence.plugins.confluence-view-file-macro:view-file-macro-editor-content-resources.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-view-file-macro:view-file-macro-editor-content-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/d41d8cd98f00b204e9800998ecf8427e-CDN/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/3.1.4/\_/download/batch/com.atlassian.confluence.plugins.confluence-software-blueprints:common-resources/com.atlassian.confluence.plugins.confluence-software-blueprints:common-resources.css" data-wrm-key="com.atlassian.confluence.plugins.confluence-software-blueprints:common-resources" data-wrm-batch-type="resource" media="all"> <link type="text/css" rel="stylesheet" href="/s/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/22/\_/styles/custom.css" media="all">

search

recentlyviewed

attachments

weblink

advanced

image-effects

image-attributes

*   [Paragraph](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   *   [Paragraph](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Heading 1](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Heading 2](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Heading 3](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Heading 4](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Heading 5](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Heading 6](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Preformatted](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Quote](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)


*   [Bold](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Italic](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Underline](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Colour picker](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   [More colours](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Formatting](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   *   [Strikethrough](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Subscript](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Superscript](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Monospace](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   *   [Clear formatting](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)


*   [Bullet list](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Numbered list](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Task list](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Outdent](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Indent](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Align left](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Align center](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Align right](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Page layout](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Link](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Table](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)


*   [Insert](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   Insert content

   *   [Files and images](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Link](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Symbol](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Emoticon](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Markup](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Horizontal rule](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Task list](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   Insert macro

   *   [User mention](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [JIRA Issue/Filter](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Info](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Add Lucidchart Diagram](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Gliffy Diagram](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Status](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Gallery](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Table of Contents](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Team Calendar](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Other macros](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)


*   [Page layout](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

   *   [No layout](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Two column (simple)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Two column (simple, left sidebar)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Two column (simple, right sidebar)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Three column (simple)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Two column](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Two column (left sidebar)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Two column (right sidebar)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Three column](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
   *   [Three column (left and right sidebars)](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)


*   [Undo](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)
*   [Redo](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Find/Replace](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

*   [Keyboard shortcuts help](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4#)

Edit

Preview

Save

Cancel

Discard unpublished changes
---------------------------

This will restore the editor to the last published version of this page. If anyone else has unpublished changes on the page, you'll be discarding those changes too.

Key: This line was added. This line was removed. Formatting was changed.

* * *

Discard See changes Close

Discard Close

<meta name="ajs-use-watch" content="true"> <meta name="ajs-attachment-source-content-id" content="26187384"> <meta name="ajs-fallback-mode" content="false"> <meta name="ajs-use-inline-tasks" content="true"> <meta name="ajs-heartbeat" content="true"> <meta name="ajs-action-locale" content="en\_GB"> <meta name="ajs-editor-plugin-resource-prefix" content="/s/en\_GB/6441/7d28db648b08c18ff79fe305b2da56665a84a1fc/5.10.4/_"> <meta name="ajs-existing-draft-id" content="0"> <meta name="ajs-user-watching-own-content" content="true"> <meta name="ajs-new-page" content="false"> <meta name="ajs-content-id" content="26187384"> <meta name="ajs-editor-mode" content="richtext"> <meta name="ajs-form-name" content="inlinecommentform"> <meta name="ajs-auto-start" content="false"> <meta name="ajs-can-attach-files" content="false"> <meta name="ajs-show-draft-message" content="false"> <meta name="ajs-conf-revision" content="confluence$content$26187384.0"> <meta name="ajs-draft-id" content="0"> <meta name="ajs-shared-drafts" content="false"> <meta name="ajs-min-editor-height" content="150"> <meta name="ajs-draft-share-id" content=""> <meta name="ajs-content-type" content="page"> <meta name="ajs-version-comment" content=""> <meta name="ajs-draft-type" content="page"> <meta name="ajs-synchrony-token" content=""> <meta name="ajs-synchrony-service-uri" content=""> <meta name="ajs-synchrony-resources-uri" content=""> <meta name="ajs-synchrony-app-id" content=""> <meta name="ajs-synchrony-dark-enabled" content="false"> <meta name="ajs-synchrony-expiry" content=""> <meta name="ajs-max-thumb-width" content="300"> <meta name="ajs-max-thumb-height" content="300"> <meta name="ajs-can-send-email" content="true"> <meta name="ajs-is-dev-mode" content="false"> <meta name="ajs-draft-save-interval" content="30000"> <meta name="ajs-show-hidden-user-macros" content="false"> <meta name="ajs-is-admin" content="false"> <meta name="ajs-confluence.prefs.editor.disable.autocomplete" content="false"> <meta name="ajs-confluence.prefs.editor.disable.autoformat" content="false"> <meta name="ajs-heartbeat-interval" content="30000"> <form id="tinymce-table-form" class="aui"> <div class="field-group"> <label for="rows">Rows</label> <input id="rows" name="rows" type="text" size="3" autocomplete="off" value="{0}"> </div> <div class="field-group"> <label for="cols">Columns</label> <input id="cols" name="cols" type="text" size="3" autocomplete="off" value="{1}"> </div> <div class="field-group hidden"> <input id="width" type="hidden" name="width" value=""> <label for="width">Width</label> </div> <div class="group"> <div class="checkbox"> <input id="table-heading-checkbox" class="checkbox" type="checkbox" name="heading" checked="checked" value="true"> <label for="table-heading-checkbox">First row is heading</label> </div> </div> <div class="group hidden"> <div class="checkbox"> <input id="table-equal-width-columns-checkbox" class="checkbox" type="checkbox" name="equal-width-columns" value="false"> <label for="table-equal-width-columns-checkbox">Equal width columns</label> </div> </div> </form>





## Overview ##
The Overview section should provide information on how the term is used in the GDC. Subsections can be created for this section, if applicable.
## References ##
1. [reference name] (reference URL)

## External Links ##
* [external link name] (external link URL)

Categories: General
