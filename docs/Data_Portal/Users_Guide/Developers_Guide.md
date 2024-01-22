# GDC Portal V2 Application Developer Guide

## Introduction

This guide will detail the process of developing applications for the GDC Portal Version 2.0. It describes the
structure of the GDC Portal, how to use the GDC Portal API, and how to develop applications for the GDC Portal.

The GDC Portal is designed to support the development of applications that allow for analysis, visualization,
and refinement of cohorts. The GDC Portal is built on top of the [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/),
which provides access to the GDC data. The GDC Portal provides a framework for developing applications that
can be used to analyze and visualize data from the GDC.

## Table of Contents

- [Introduction](#introduction)
  - [Overview of an Application](#overview-of-an-application)
  - [Local vs Global Filters](#local-vs-global-filters)
  - [Cohorts and Filters](#cohorts-and-filters)
- [Using the Portal Application API](#using-the-portal-application-api)
  - [Case Information](#case-information)
  - [File Information](#file-information)
  - [Sets: Gene, SSMS, and Case](#sets-gene-ssms-and-case)
  - [Creating a cohort](#creating-a-cohort)
  - [Altering a cohort](#altering-a-cohort)
    - [Updating, removing, and clearing filters](#updating-removing-and-clearing-filters)
    - [Updating the cohort name](#updating-the-cohort-name)
    - [Setting the current cohort](#setting-the-current-cohort)
  - [Count Information](#count-information)
  - [Component Library](#component-library)
    - [Buttons](#buttons)
    - [Modals](#modals)
    - [Charts](#charts)
    - [Facets](#facets)
    - [VerticalTable](#verticaltable)
- [Application Development](#application-development)
  - [Getting Started](#getting-started)
  - [Application Layout](#application-layout)
  - [Local State](#local-state)
  - [Persisting the local State](#persisting-the-local-state)
  - [Application Hooks](#application-hooks)
  - [Creating a new cohort](#creating-a-new-cohort)
  - [Application Registration](#application-registration)
  - [Source code layout](#source-code-layout)


- [Appendix](#appendix)
  - [Using selectors and hooks](#using-selectors-and-hooks)
    - [Selectors](#selectors)
    - [Hooks](#hooks)
  - [Querying the GDC API Directly](#querying-the-gdc-api-directly)

## Introduction

The GDC Portal is designed to support the development of applications that allow for analysis, visualization,
and refinement of cohorts. The GDC Portal is built on top of the [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/),
which provides access to the GDC data. The GDC Portal provides a framework for developing applications that
can be used to analyze and visualize data from the GDC. The GDC Portal is built on top of the [React](https://reactjs.org/) framework and uses the [Redux](https://redux.js.org/) library for state management. The GDC Portal uses
NextJS to provide server-side rendering of React components. Mantine.dev is the base component library used
and styling is done with [TailwindCSS](https://tailwindcss.com/).

The GDC Portal contains an Analysis Center where applications are displayed for users to use with their cohorts.
The GDC Portal also provides a framework for developing applications that can be used to analyze and visualize data from the GDC.

## Overview of an Application

Applications are High Order Components (HOC) that are rendered in the Analysis Center. The portal major functions
like Project, Downloads, and Protein Paint are all applications. Each application handles a specific task and can be used to
refine and analyze cohorts. Applications have access to all the current cohort information and can use that information
to query the GDC API for additional information.

Local and Cohort filters are available to applications. Local filters are filters that are specific to the application and are used to refine the data that is displayed in the application. Local filters are those available from the GDC API and are typically not the
most common. For example in the Mutation Frequency application, the local filters are the gene and mutation type filters. In the figure
below the local filters are highlighted in yellow. These filters are used to refine the input cohort allowing users to
drill down to specific genes and mutation types of interest in the cohort.

![Mutation Frequency](./images/mutation_frequency_app.png)

### Local vs Global Filters

A Portal application's input can be anything including a single cohort or multiple cohorts. The application then can either add filters to refine the cohort by adding filters, create additional cohorts, or display the data in a visualization. Applications typically
have:
* local filters that are used to refine the data displayed in the application.
* UI components that are used to display the data in the application.
* State that is used to store the data displayed in the application.
* Actions that are used to update the state of the application.

Applications can also create new cohorts. These cohorts can be used by other Portal applications.

![Structure of an Application](./images/application_structure.png)

## Cohorts and Filters

From an application perspective, a cohort is an Object containing the following information:
```typescript
interface Cohort {
  id: string;        // unique id for cohort
  name: string;      // name of cohort
  filters: FilterSet; // active filters for cohort
  caseSet: CaseSetDataAndStatus; // case ids for frozen cohorts
  modified?: boolean; // flag which is set to true is modified and unsaved
  modified_datetime: string; // last time cohort was modified
  saved?: boolean; // flag indicating if cohort has been saved.
  counts: CountsDataAndStatus; //case, file, etc. counts of a cohort
}
```

Likely the most important part of the cohort is the `filters` field. The `filters` field contains the active filters for the cohort.
The `filters` field is a `FilterSet` object. The `FilterSet` object contains the active filters for the cohort. When calling either the
GDC API or GDC GraphQL API the `FilterSet`` is converted to the appropriate format for the API. The `FilterSet` object is of the form:

```typescript
interface FilterSet {
  op: "and" | "or"; // operator for combining filters
  root: Record<string,Operation >; // map of filter name to filter operation
}
```
Operation are GDC API filters described in the [GDC API Guide](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query). These are:
* Equals
* NotEquals
* LessThan
* LessThanOrEquals
* GreaterThan
* GreaterThanOrEquals
* Exists
* Missing
* Includes
* Excludes
* ExcludeIfAny
* Intersection
* Union

The `root` field is a map of filter names (as defined in the GDC API) to filter operation. The filter operation can be either a single operation
or a `FilterSet` object. The `op` field will eventually support either `and` or `or`, however at this time only `and` is supported. The `and`operator is used to combine filters using the `and` operator. The `or` operator is used to combine
filters using the `or` operator. The `FilterSet` object is converted to the appropriate format for the GDC API when the cohort is saved.
When using the GDC REST API, the FilterSet can be converted into the appropriate format using the `filterSetToOperation` function.
When using the GDC GraphQL API, the FilterSet can be using the `convertFilterSetToGraphQL` function. The API guide will provide information
on what format the filters should be in for the API. Also as the code is in TypeScript, the IDE will provide information on the format as well.


### Getting Cohort Information

The current active cohort can be accessed via the selector `selectCurrentCohort`. This selector returns the current cohort,
which is the cohort that is currently being displayed in the Cohort Management Bar. Accessing the current cohort is done via the
selector:

```typescript
import {useCoreSelector,  selectCurrentCohort } from '@gff/core';

const currentCohort = useSelector(selectCurrentCohort);
```

By using the selector, the component/application will be updated when the cohort changes. There are also selectors for getting a particular field from the cohort. For example, to get the cohort name, the selector `selectCurrentCohortName` can be used. The selectors are:


* `selectCurrentCohort`
* `selectCurrentCohortName`
* `selectCurrentCohortId`
* `selectCurrentCohortFilters`
* `selectCurrentCohortModified`
* `selectCurrentCohortModifiedDatetime`
* `selectCurrentCohortSaved`
* `selectCurrentCohortCounts`

The current active filters can be accessed via the selector `selectCurrentCohortFilters`. This selector returns the current filters,
which are the filters that are currently being displayed in the Cohort Management Bar. Accessing the current filters is done via the
selector:

```typescript
import {useCoreSelector,  selectCurrentFilters } from '@gff/core';

const currentFilters = useSelector(selectCurrentCohortFilters);
```

By using the selector, the application will be updated when the filters change. The filters are returned as a `FilterSet` object described above.

All the cohorts can be selected using the selector `selectAllCohorts`. This selector returns all the cohorts in the store. Accessing all the cohorts is done via the selector:

```typescript
import {useCoreSelector,  selectAllCohorts } from '@gff/core';

const allCohorts = useSelector(selectAllCohorts);
```
# Using the Portal Application API

The GDC Portal provides a number of hooks for querying the GDC API. These hooks are located in the `@gff/core` package.
The hooks are designed to work in a manner similar to the RTL Query hooks. The hooks take arguments and return an object.
The object contains the data and the status of the query. The status of the query is stored in the `isSuccess` variable.
The @gff/core package also provides a set of selectors that return values stored in the core redux store: `CoreStore`.

There are a number of hooks and selectors that are available for querying the GDC API, a subset of which are shown below:

![hooks and selectors](./images/hooks_and_selectors.png)

## Case Information

The GDC Portal provides several hooks for querying case information. These hooks are located in the `@gff/core` package.
Cases can be queried using several different methods. The 'useAllCases' hook returns all the cases in the GDC and can be filtered by the current cohort as shown below:

```typescript
import { useCoreSelector, useAllCases } from '@gff/core';

...

const [pageSize, setPageSize] = useState(10);
const [offset, setOffset] = useState(0);
const [searchTerm, setSearchTerm] = useState<string>("");
const [sortBy, setSortBy] = useState<SortBy[]>([]);
const cohortFilters = useCoreSelector((state) =>
        selectCurrentCohortFilters(state),
);


const { data, isFetching, isSuccess, isError, pagination } = useAllCases({
  fields: [
    "case_id",
    "submitter_id",
    "primary_site",
    "disease_type",
    "project.project_id",
    "project.program.name",
    "demographic.gender",
    "demographic.race",
    "demographic.ethnicity",
    "demographic.days_to_death",
    "demographic.vital_status",
    "diagnoses.primary_diagnosis",
    "diagnoses.age_at_diagnosis",
    "summary.file_count",
    "summary.data_categories.data_category",
    "summary.data_categories.file_count",
    "summary.experimental_strategies.experimental_strategy",
    "summary.experimental_strategies.file_count",
    "files.file_id",
    "files.access",
    "files.acl",
    "files.file_name",
    "files.file_size",
    "files.state",
    "files.data_type",
  ],
  size: pageSize,
  filters: cohortFilters,
  from: offset * pageSize,
  sortBy: sortBy,
  searchTerm,
});

```
The `useAllCases` hook takes a number of arguments:
* `fields` - the fields to return from the GDC API
* `size` - the number of cases to return
* `filters` - the filters to apply to the cases
* `from` - the starting index of the cases to return
* `sortBy` - the fields to sort the cases by
* `searchTerm` - the search term to use to search the cases

This call is used in the Table view tab of the Cohort Management Bar.

Information for a single case can be queried using the `useCaseSummary` hook. This call is used in the caseView page:
[portal.gdc.cancer.gov/cases/5693302a-4548-4c0b-8725-0cb7c67bc4f8](https://portal.gdc.cancer.gov/cases/5693302a-4548-4c0b-8725-0cb7c67bc4f8)


```typescript

  const { data, isFetching } = useCaseSummary({
  filters: {
    content: {
      field: "case_id",
      value: case_id,
    },
    op: "=",
  },
  fields: [
    "files.access",
    "files.acl",
    "files.data_type",
    "files.file_name",
    "files.file_size",
    "files.file_id",
    "files.data_format",
    "files.state",
    "files.created_datetime",
    "files.updated_datetime",
    "files.submitter_id",
    "files.data_category",
    "files.type",
    "files.md5sum",
    "case_id",
    "submitter_id",
    "project.name",
    "disease_type",
    "project.project_id",
    "primary_site",
    "project.program.name",
    "summary.file_count",
    "summary.data_categories.file_count",
    "summary.data_categories.data_category",
    "summary.experimental_strategies.experimental_strategy",
    "summary.experimental_strategies.file_count",
    "demographic.ethnicity",
    "demographic.demographic_id",
    "demographic.gender",
    "demographic.race",
    "demographic.submitter_id",
    "demographic.days_to_birth",
    "demographic.days_to_death",
    "demographic.vital_status",
    "diagnoses.submitter_id",
    "diagnoses.diagnosis_id",
    "diagnoses.classification_of_tumor",
    "diagnoses.age_at_diagnosis",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.days_to_last_known_disease_status",
    "diagnoses.days_to_recurrence",
    "diagnoses.last_known_disease_status"]
  });
```
The `useCaseSummary` hook takes a number of arguments:
* `fields` - the fields to return from the GDC API
* `filters` - the filters to apply to the cases and where the caseId is passed in

## File Information

Similar to the case information, the GDC Portal provides a number of hooks for querying file information. These hooks are located in the `@gff/core` package.
To get a list of files associated with a cohort, the `useGetFilesQuery` hook can be used. This call is used in the repository application. The repository application is used to display the files associated with a cohort. The `useGetFilesQuery` hook takes a number of arguments:

```typescript
import {
  useCoreDispatch,
  useCoreSelector,
  selectCurrentCohortFilters,
  buildCohortGqlOperator,
  joinFilters,
  useFilesSize,
} from "@gff/core";

...

const coreDispatch = useCoreDispatch();
const [sortBy, setSortBy] = useState<SortBy[]>([]); // states to handle table sorting and pagination
const [pageSize, setPageSize] = useState(20);
const [offset, setOffset] = useState(0);

const repositoryFilters = useAppSelector((state) => selectFilters(state)); // as this is a app get the repository filters from the app state (local filters)
const cohortFilters = useCoreSelector((state) =>    // get the cohort filters from the core state (global filters)
        selectCurrentCohortFilters(state),
);

const { data, isFetching, isError, isSuccess } = useGetFilesQuery({
  case_filters: buildCohortGqlOperator(cohortFilters),
  filters: buildCohortGqlOperator(repositoryFilters),
  expand: [
    "annotations", //annotations
    "cases.project", //project_id
    "cases",
  ],
  size: pageSize,
  from: offset * pageSize,
  sortBy: sortBy,
});

```

The `useGetFilesQuery` hook takes a number of arguments:
* `case_filters` - the filters to apply to the cases
* `filters` - the filters to apply to the files
* `expand` - the fields to expand
* `size` - the number of files to return
* `from` - the starting index of the files to return
* `sortBy` - the fields to sort the files by

Note this hook was designed to take global filters (e.x the current cohort as `case_filters`) and local filters (the repository filters).

Information for a single file can be queried using the `useFileSummary` hook. This call is used in the `File Summary View` page [portal.gdc.cancer.gov/files/0b5a9e7e-8e2e-4b7a-9b7e-ff5d9c5b2b2b](https://portal.gdc.cancer.gov/files/0b5a9e7e-8e2e-4b7a-9b7e-ff5d9c5b2b2b)

```typescript
 const { data: { files } = {}, isFetching } = useGetFilesQuery({
  filters: {
    op: "=",
    content: {
      field: "file_id",
      value: setCurrentFile,
    },
  },
  expand: [
    "cases",
    "cases.annotations",
    "cases.project",
    "cases.samples",
    "cases.samples.portions",
    "cases.samples.portions.analytes",
    "cases.samples.portions.slides",
    "cases.samples.portions.analytes.aliquots",
    "associated_entities",
    "analysis",
    "analysis.input_files",
    "analysis.metadata.read_groups",
    "downstream_analyses",
    "downstream_analyses.output_files",
    "index_files",
  ],
});
```
The `useFileSummary` hook takes several arguments:
* `filters` - the filters to apply to the cases and where the file uuid is passed in
* `expand` - the fields to expand

## Sets: Gene, SSMS, and Case

Sets are supported by the GDC API and are used to create an entity that represents a set of items as a `set_id`. Sets can be
gene sets, SSM sets, or case sets. All GDC APIs support passing sets as a filter parameter.
The GDC Portal provides a number of hooks for creating and querying set information.

a set can be created using one of the following hooks:

* `useCreateGeneSetFromValuesMutation`
* `useCreateSsmsSetFromValuesMutation`
* `useCreateCaseSetFromValuesMutation`
* `useCreateGeneSetFromFiltersMutation`
* `useCreateSsmsSetFromFiltersMutation`
* `useCreateCaseSetFromFiltersMutation`

These functions will create a set from either a list of values or a filter set. The create from Values hooks take a single
parameter `values` which is an array of values, while the create from filters hooks take one required parameter `filters`
which is a filter set or JSON object. Both calls return the created `set_id` if the set was successfully created.

As the above hooks are Redux Toolkit Query hooks, namely mutation hooks, they return a tuple of the form:
`[mutationHook, response]` which is a function to call the mutation and the response from the mutation. The mutation hook can be used like:

```typescript
 const [createSet, response] = createSetHook();

  const handleCreateSet = async () => {
    const { data } = await createSet({
      variables: {
        values: ["TP53", "KRAS", "EGFR"],
      },
    });
    if (response.isSuccess) {
      dispatch(
        addSet({
          setType,
          setName: form.values.name.trim(),
          setId: response.data as string,
        }),
      );
    }
    ;
  }
 ```

Once a set is created it can be altered using the following hooks:
* `useAppendToGeneSetMutation`
* `useAppendToSsmSetMutation`
* `useRemoveFromGeneSetMutation`
* `useRemoveFromSsmSetMutation`

Sets can be managed using the following actions:

* `addSet`
* `removeSet`
* `updateSet`

The following selectors are available for getting set information:

* `selectAllSets`
* `selectSetById`
* `selectSetByName`
* `selectSetByType`

Finally, the following hooks are available for querying set size:

* `useGeneSetCountsQuery`
* `useSsmSetCountsQuery`
* `useCaseSetCountsQuery`

## Creating a cohort

Depending on what your application does, you may want to create a new cohort. Although the GDC Portal SDK provides a
number of functions for creating a new cohort. It is highly recommended that the application use the provided Button and
SaveCohortModal components to create a new cohort. The Button and SaveCohortModal components are located in
the `@gff/portal-proto` package.

To create a cohort using the SaveCohortModal component the following code can be used:
In summary, the above code flow is:

1. The ProjectsCohortButton component renders a button with the label "Save New Cohort".
2. When the button is clicked, it sets the state variable `showSaveCohort` to true, which triggers the rendering of the `SaveCohortModal` component.
3. The `SaveCohortModal` component passed:
   * an onClose function that sets the showSaveCohort state variable to false.
   * a filters prop, which is an object defining the filters for the cohort based on the selected projects.
4. The `SaveCohortModal` will use the passed filter to create, name and save the cohort when the save button is clicked.

Additional details on the `SaveCohortModal` component can be found in the [Component Library](#component-library) section as
well as buttons to create a saved cohort.


## Altering a cohort

Altering a cohort is done by dispatching actions o add, remove, or clear filters. The following actions are available
for altering the current cohort:

* `updateCohortFilter`
* `removeCohortFilter`
* `clearCohortFilters`

Note that all of these operations are applied to the current cohort. The current cohort is the cohort that is currently
being displayed in the Cohort Management Bar. The current cohort can be accessed via the `selectCurrentCohort` selector.
The current cohort's filters can be accessed via the `selectCurrentCohortFilters` selector.

### Updating, removing, and clearing filters

to update the current selected cohort's filter, the `updateCohortFilter` action can be used. The `updateCohortFilter` action takes two arguments:
```typescript

interface UpdateFilterParams {
  field: string;
  operation: Operation;
}

```
where `field` is the field to update and `operation` is the operation to apply to the field. For example to update the
`cases.project.project_id` field to include the project `TCGA-ACC` the following code can be used:

```typescript
import { useCoreDispatch, updateCohortFilter } from '@gff/core';

const coreDispatch = useCoreDispatch();

coreDispatch(updateCohortFilter({
  field: "cases.project.project_id",
  operation: {
    op: "in",
    content: {
      field: "cases.project.project_id",
      value: ["TCGA-ACC"],
    },
  },
}));
```

This will update the current cohort's filter to include the project `TCGA-ACC`. The `removeCohortFilter` action can be used to remove a filter from the current cohort. The `removeCohortFilter` action takes a single argument:
```typescript

interface RemoveFilterParams {
  field: string;
}

```
where `field` is the field to remove. For example to remove the `cases.project.project_id` field from the current cohort's filter the following code can be used:

```typescript
import { useCoreDispatch, removeCohortFilter } from '@gff/core';

const coreDispatch = useCoreDispatch();

coreDispatch(removeCohortFilter({
  field: "cases.project.project_id",
}));
```

This will remove the `cases.project.project_id` field from the current cohort's filter. The `clearCohortFilters` action can be used to clear all the filters from the current cohort. The `clearCohortFilters` action takes no arguments. For example to clear all the filters from the current cohort the following code can be used:

```typescript
import { useCoreDispatch, clearCohortFilters } from '@gff/core';

const coreDispatch = useCoreDispatch();

coreDispatch(clearCohortFilters());
```

This will clear all the filters from the current cohort.

#### Updating the cohort name

The cohort name can be updated using the `updateCohortName` action. The `updateCohortName` action takes a single argument:
```typescript

interface UpdateCohortNameParams {
  name: string;
}

```
where `name` is the new name for the cohort. For example, to update the current cohort's name to `My Cohort` the following code can be used:

```typescript
import { useCoreDispatch, updateCohortName } from '@gff/core';

const coreDispatch = useCoreDispatch();

coreDispatch(updateCohortName({
  name: "My Cohort",
}));
```

This will update the current cohort's name to `My Cohort`.

#### Setting the current cohort

The current cohort can be set using the `setCurrentCohort` action. The `setCurrentCohort` action takes a single argument:
```typescript

interface SetCurrentCohortParams {
  cohortId: string;
}

```
where `cohortId` is the id of the cohort to set as the current cohort. For example to set the cohort with id `1234` as the current cohort the following code can be used:

```typescript
import { useCoreDispatch, setCurrentCohort } from '@gff/core';

const coreDispatch = useCoreDispatch();

coreDispatch(setCurrentCohort({
  cohortId: "1234",
}));
```

This will set the cohort with ID `1234` as the current cohort.


## Count Information

Counts information can be queried using the `useTotalCounts` hook. This hook takes a number of arguments:

```typescript
import { useTotalCounts } from "@gff/core";

const { data, isFetching, isSuccess, isError } = useTotalCounts();
```

this will return the total counts for the GDC. The data in the response is of the form:

```typescript
interface TotalCounts {
  counts: {
    caseCounts: number;
    fileCounts: number;
    genesCounts: number;
    mutationCounts: number;
    repositoryCaseCounts: number;
    projectsCounts: number;
    primarySiteCounts: number;
  },
  status: DataStatus;
}
```
where `DataStatus` is defined as:

```typescript
export type DataStatus = "uninitialized" | "pending" | "fulfilled" | "rejected";
```

## Component Library

As a developer you will likely want to use the components provided by the GDC Portal. The GDC Portal provides a number of components
that make it easy to develop applications for the GDC Portal. These components are located in the `@gff/portal-proto` package.
In several components, the GDC Portal uses the [Mantine](https://mantine.dev/) component library but base components and encapsulates calls to the GDC API so that you do not have to.

### Buttons

The GDC Portal provides a number of buttons that can be used for various purposes. These buttons are located in the `@gff/portal-proto` package.
The buttons are:

* `DownloadButton` - a button that can be used to download data from the GDC API.
* `SaveCohortButton` - a button that can be used to save a cohort.

The `DownloadButton` component is used in the repository application to download data from the GDC API. The `DownloadButton` component takes a number of arguments:

```tsx

<DownloadButton
        inactiveText={`Download ${numFilesCanAccess} Authorized File${
                numFilesCanAccess !== 1 ? "s" : ""
        }`}
        activeText=""
        disabled={
                numFilesCanAccess === 0 ||
                (user.username && dbGapList.length > 0 && !checked)
        }
        endpoint="data"
        extraParams={{
          ids: (filesByCanAccess?.true || []).map((file) => file.file_id),
          annotations: true,
          related_files: true,
        }}
        method="POST"
        setActive={setActive}
/>

```

The parameters for the DownloadButton are defined in the Portal V2 SDK API documentation. The `DownloadButton` component will take care of
calling the GDC API and downloading the data. The `DownloadButton` component will also provide status that can be used with
a progress bar or spinner to display the progress of the download.

The `SaveCohortButton` component is used in the repository application to save a cohort.
![img.png](images/create_cohort_button.png)
The `SaveCohortButton` component takes a number of arguments:

```tsx

<CohortCreationButton
        numCases={cohort1Count}
        label={cohort1Count.toLocaleString()}
        filtersCallback={async () =>
                generateFilters(caseSetIds[0], caseSetIds[1])
        }
/>
```

The CohortCreationButton component will show the number of selected cases and will create a new saved cohort
when the button is clicked. The `filtersCallback` is a function that returns the filters for the cohort.

### Modals

Modals are used to show transitory information or obtain information from the user. The GDC Portal provides many
modals that can be used for various purposes. One such modal is the `SaveCohortModal` component mentioned previously.

![img.png](images/save_cohort_modal.png)

* `SaveCohortModal` - a modal that can be used to save a cohort.
* Various modals for displaying information on Sets:
  * `CaseSetModal`
  * `GeneSetModal`
  * `MutationSetModal`
* `SaveOrCreateEntityModal` - a modal that can be used to save or create a new entity.

These modals and others, are documented in the Portal V2 SDK API documentation.

### Charts

Basic charts are provided for use in your application, although applications are free to use any charting library they wish.
The charts provided are:

* `BarChart` - a bar chart

  ![Bar Chart](images/primary_site.png)

The `BarChart` component (based on Plotly) is passed data in the form:

```typescript
import { PlotData } from "plotly.js";


export interface BarChartData {
  datasets: Partial<PlotData>[];
  yAxisTitle?: string;
  tickvals?: number[];
  ticktext?: string[];
  label_text?: string[] | number[];
  title?: string;
  filename?: string;
}
```

Where `datasets` is an array of PlotData objects, which at the minimum contain the `x` and `y` fields.

The `yAxisTitle` is the title for the y-axis, the `tickvals` and `ticktext` are the tick values and text for
the x-axis, the `label_text` is the text for the labels, the `title` is the title for the chart,
and the `filename` is the filename to use when downloading the chart.

Note that BarChart needs to be imported as a dynamic component:

```tsx
import dynamic from "next/dynamic";

const BarChart = dynamic(() => import("@/components/charts/BarChart"), {
  ssr: false,
});
```

* `Cancer Distribution` - a cancer distribution chart

  ![cancer distribution](images/most-frequently-mutated-genes-bar-chart.png)

The `CancerDistribution` component (based on Plotly) is different as it passed the Gene Symbol
and optionally cohort and gene filters.

```typescript
interface CNVPlotProps {
  readonly gene: string;
  readonly height?: number;
  readonly genomicFilters?: FilterSet;
  readonly cohortFilters?: FilterSet;
}
```

These charts and others are documented in the Portal V2 SDK API documentation.

### Facets

Facet components are provided for use in building local filters for your application. There are two types of facet components:

* `EnumFacet` - a facet that is used to filter on an enum field
* `DateFacet` - a facet that is used to filter a date field
* `NumericRangeFacet` - a facet that is used to filter on a range field
* `PercentileFacet` - a facet that is used to filter on a percentile field
* 'AgeRangeFacet' - a facet that is used to filter on an age range field
* `TextFacet` - a facet that is used to filter a text field
* `BooleanFacet` - a facet that is used to filter on a boolean field

<img src="images/components/enum_facet.png" alt="Enum Facet Component" width="800" height="auto">

*Enum Facet*

<img src="images/components/numeric_range_facet.png" alt="Range Facet Component" width="600" height="auto">

*Range Facet*

<img src="images/components/date_range_facet.png" alt="Date Range Facet Component" width="600" height="auto">

*Date Range Facet*


<img src="images/components/number_range.png" alt="Number Range Facet Component" width="600" height="auto">

*Number Range Facet*

<img src="images/components/percentile_facet.png" alt="Percent Range Facet Component" width="600" height="auto">

*Percentile Facet*


<img src="images/components/age_range_facet.png" alt="Age Range Facet Component" width="600" height="auto">

*Age Range Facet*


<img src="images/components/exact_value_facet.png" alt="Exact Value Facet Component" width="600" height="auto">

*Exact Value Facet*

<img src="images/components/toggle_facet.png" alt="Boolean Toggle Facet Component" width="600" height="auto">

*Toggle Facet*


The facet components are documented in the Portal V2 SDK API documentation. As these components are passed data fetcher
and filter management hooks, they can be used for both cohort and local filters in an application.

### VerticalTable

The VerticalTable component is used to display data in a table format. The VerticalTable component is a Mantine component
implementing react table version 8. The VerticalTable component has a number of parameters, the most important
being data, columns, and filters. The data is the data to display in the table, the columns are the columns to display
in the table, and is where the fields of the data tp be rendered are.
The table has support for searching, sorting, and pagination. It can be configured to render many different types of columns, including text, numeric, and date. The table can also be configured
to use React components for rendering columns.
The Vertical Table is used for most of the table views in the GDC Portal. There are a number of examples of its use and
is documented in the Portal V2 SDK API documentation.

![vertical_table.png](images%2Fcomponents%2Fvertical_table.png)
*Vertical Table*

# Application Development

## Getting Started

The GDC Portal V2 is a monorepo that contains all the code for the GDC Portal. The monorepo is managed using [lerna](https://lerna.js.org) and [npm]().
The monorepo contains the following packages:

* `@gff/core` - contains the core components and hooks for the GDC Portal.
* `@gff/portal-proto` - contains the UI components and application framework (using NextJS) for the GDC Portal.

Note that the UI components located in the `@gff/portal-proto` package will be refactored into a separate package in the future, and
`@gff/portal-proto` will be renamed to `@gff/portal`.

You can get started by cloning the repo and following the instructions in the [README.md](https://github.com/NCI-GDC/gdc-frontend-framework/blob/develop/README.md) file.

## Application Layout

A typical application will have the following layout. The main section of the application is the area where components like tables, graphs, and other components are displayed. Local filters are displayed on the left side and
depending on the numbers will scroll vertically. This is a typical layout but other layouts are possible, like
in the case of Protein Paint. Applications are encouraged to use vertical space as much as possible, as horizontal
scrolling can be a poor user experience.

This section will describe parts of the Project application and how it is structured. The Project application is a
simple application that displays a table of projects and allows the user to filter the projects by a number of filters.
As the local filters are selected the table display is updated, but the cohort is not changed (e.i cohort filters are not
updated). The Project application is located in the `@gff/portal-proto` package in the `src/features/projectsCenter` directory.
The user can create a new saved cohort by selecting projects and clicking the "Save New Cohort" button. This will open
a modal that will allow the user to name the cohort and save it. The Project application is a good example of how to use
the GDC Portal SDK to create an application.

![projects_app_parts_outlined.png](images%2Fprojects_app_parts_outlined.png)
*Major Sections of an Application*

## Local State

Depending on the application, it may be necessary to maintain the local state. For example, in the Projects application,
the selected local filters, in this case, represented as Enumeration Facets, are stored in the local state. This allows
the application to remember the selected filters when the user navigates away from the page and then returns. Persisting
the state uses [Redux Toolkit] and [Redux Persist] to store the state in local storage. While the CoreState is managed
by the portal core, the local state is managed by the application. Using a separate store for the local state allows
the application to manage the state without having to worry about affecting the core state.

The Portal core provides a number of functions to assist in the creation and persisting of the redux store and will create
handlers such as AppState, AppDispatch, and AppSelector. The AppState is the type of the local state, the AppDispatch
is the type of the dispatch function, and the AppSelector is the type of the selector function.

An application can create all of the them using the `createAppStore` function:

```typescript
import { createAppStore } from "@gff/core";
import { projectCenterFiltersReducer } from "./projectCenterFiltersSlice";

const PROJECT_APP_NAME = "ProjectCenter";

// create the store, context and selector for the ProjectsCenter
// Note the project app has a local store and context which isolates
// the filters and other store/cache values

const reducers = combineReducers({  
  projectApp: projectCenterFiltersReducer, // Your application might have more that one reducer
});

export const { id, AppStore, AppContext, useAppSelector, useAppDispatch } =
        createAppStore({
          reducers: reducers,
          name: PROJECT_APP_NAME,
          version: "0.0.1",
        });

export type AppState = ReturnType<typeof reducers>;
```
Call this function will create the local store given the reducers and the name, and version of the application.
The name of the application is used to create the local storage key for the application. The `id` is the id of the
application and is used to create the local storage key. The `AppStore` is the local store, the `AppContext` is the
local context, the `useAppSelector` is the selector hook, and the `useAppDispatch` is the dispatch hook.

Since you now have a local store, you can create a slice for the local state. The slice is a standard Redux Toolkit
slice and will contain the reducer, actions, and selectors for the local state.

## Persisting the local state

If is desirable to persist the local state. This can be done using the `persistReducer` function from the
[redux-persist](https://github.com/rt2zz/redux-persist) package. Any reducer can be persisted by creating a
persisted store and passing the reducer to the `persistReducer` function. For example, the `createAppStore` function
can be modified to persist the local filter state as:



```typescript #title="src/appApi.ts"
import { combineReducers } from "redux";
import { persistReducer } from "redux-persist";
import storage from "redux-persist/lib/storage";
import { createAppStore } from "@gff/core";
import { projectCenterFiltersReducer } from "./projectCenterFiltersSlice";

const PROJECT_APP_NAME = "ProjectCenter";

const persistConfig = {
  key: PROJECT_APP_NAME,
  version: 1,
  storage,
  whitelist: ["projectApp"],
};

// create the store, context and selector for the ProjectsCenter
// Note the project app has a local store and context which isolates
// the filters and other store/cache values

const reducers = combineReducers({
  projectApp: projectCenterFiltersReducer,
});

export const { id, AppStore, AppContext, useAppSelector, useAppDispatch } =
  createAppStore({
    reducers: persistReducer(persistConfig, reducers),
    name: PROJECT_APP_NAME,
    version: "0.0.1",
  });

export type AppState = ReturnType<typeof reducers>;

```

For example, the `projectCenterFiltersSlice.ts` which handles the local filters, is defined as:

```typescript
import { createSlice, PayloadAction } from "@reduxjs/toolkit";
import { Operation, FilterSet } from "@gff/core";
import { AppState } from "./appApi";

export interface ProjectCenterFiltersState {
  readonly filters: FilterSet;
}

const initialState: ProjectCenterFiltersState = {
  filters: { mode: "and", root: {} },
};

const slice = createSlice({
  name: "projectCenter/filters",
  initialState,
  reducers: {
    updateProjectFilter: (
            state,
            action: PayloadAction<{ field: string; operation: Operation }>,
    ) => {
      return {
        ...state,
        filters: {
          mode: "and",
          root: {
            ...state.filters.root,
            [action.payload.field]: action.payload.operation,
          },
        },
      };
    },
    removeProjectFilter: (state, action: PayloadAction<string>) => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const { [action.payload]: _, ...updated } = state.filters.root;
      return {
        ...state,
        filters: {
          mode: "and",
          root: updated,
        },
      };
    },
    clearProjectFilters: () => {
      return { filters: { mode: "and", root: {} } };
    },
  },
  extraReducers: {},
});

export const projectCenterFiltersReducer = slice.reducer;
export const { updateProjectFilter, removeProjectFilter, clearProjectFilters } =
        slice.actions;

export const selectFilters = (state: AppState): FilterSet | undefined =>
        state.projectApp.filters;

export const selectProjectFiltersByName = (
        state: AppState,
        name: string,
): Operation | undefined => {
  return state.projectApp.filters.root[name];
};
```

The above code creates a slice for the local state. The slice contains the reducer and the actions for the local state.

The reducer is `projectCenterFiltersReducer` and the actions are `updateProjectFilter`, `removeProjectFilter`, and `clearProjectFilters`.
The selectors are `selectFilters` and `selectProjectFiltersByName`. The `selectFilters` selector returns the filters for the
application, while the `selectProjectFiltersByName` selector returns the filter for a given name.

## Application Hooks

The above can be used to define hooks for use in the local filter EnumFacet component. For example, the `useProjectFiltersByName` hook is
implemented as:

```typescript
export const useUpdateProjectsFacetFilter = (): UpdateFacetFilterFunction => {
  const dispatch = useAppDispatch();
  // update the filter for this facet

  return (field: string, operation: Operation) => {
    dispatch(updateProjectFilter({ field: field, operation: operation }));
  };
};
```

Clearing the local filters is done using the `clearProjectFilters` action:

```typescript
export const useClearProjectsFacetFilters = (): ClearFacetFiltersFunction => {
  const dispatch = useAppDispatch();
  // clear the filters for this facet

  return () => {
    dispatch(clearProjectFilters());
  };
};
```

Note that the dispatch function is the `useAppDispatch` hook returned by the `createAppStore` function. The user-selected local filters can be retrieved using the `useProjectsFilters` hook created by combining the
`useAppSelector` hook and the `selectFilters selector`:

```typescript
export const useProjectsFilters = (): FilterSet => {
  return useAppSelector((state) => selectFilters(state));
};
```
## Creating a new cohort

The project application allows users to create a new cohort from the selected projects. The cohort is created using the
`SaveCohortModal` component. The `SaveCohortModal` component passes the current cohort filters and the local project
filters to create a new saved cohort. In the case of the project application, the `SaveCohortModal` component is used
in a button component. The button component is passed the selected projects and the `SaveCohortModal` component is
rendered when the button is clicked. The `SaveCohortModal` component passes the current cohort filters and the local
project filters to create a new saved cohort. The `SaveCohortModal` component is used in the project application as:

```tsx
import React, { useState } from "react";
import { Button, Tooltip } from "@mantine/core";
import { CountsIcon } from "@/components/tailwindComponents";
import SaveCohortModal from "@/components/Modals/SaveCohortModal";

const ProjectsCohortButton = ({ pickedProjects, }: { pickedProjects: string[]; }): JSX.Element => {
  const [showSaveCohort, setShowSaveCohort] = useState(false);

  return (
          <>
            <Tooltip
                    label="Save a new cohort of cases in selected project(s)"
                    withArrow
            >
        <span>
          <Button
                  data-testid="button-create-new-cohort-projects-table"
                  variant="outline"
                  color="primary"
                  disabled={pickedProjects.length == 0}
                  leftIcon={
                    pickedProjects.length ? (
                            <CountsIcon $count={pickedProjects.length}>
                              {pickedProjects.length}{" "}
                            </CountsIcon>
                    ) : null
                  }
                  onClick={() => setShowSaveCohort(true)}
                  className="border-primary data-disabled:opacity-50 data-disabled:bg-base-max data-disabled:text-primary"
          >
            Save New Cohort
          </Button>
        </span>
            </Tooltip>
            {showSaveCohort && (
                    <SaveCohortModal
                            onClose={() => setShowSaveCohort(false)}
                            filters={{
                              mode: "and",
                              root: {
                                "cases.project.project_id": {
                                  operator: "includes",
                                  field: "cases.project.project_id",
                                  operands: pickedProjects,
                                },
                              },
                            }}
                    />
            )}
          </>
  );
};

export default ProjectsCohortButton;
```

This custom button component used the state variable `showSaveCohort` to determine if the `SaveCohortModal` component needs to be shown.
The `SaveCohortModal` component is passed the current list of projects selected by the user and handles the creation of the cohort and saving it.

### Application Demo

In addition to a application that works on cohorts, an application can have a demo. This demo can be used to show the
application's functionality. The demo is shown when the demo button is clicked. The demo button is shown when the
application is registered with `hasDemo: false,` as described in the [Application Registration](#application-registration) section.

The application can determine if the demo button should be shown by using the `useHasDemo` hook. The `useHasDemo` hook
returns a boolean indicating if the demo button should be shown. The demo button can be shown using the following code:

```tsx
import { useIsDemoApp } from "@/hooks/useIsDemoApp";

const GenesAndMutationFrequencyAnalysisTool: React.FC = () => {
  const isDemoMode = useIsDemoApp();
  ...

```

### Application Registration

An application needs to be "registered" to be used in the GDC Portal. Registration is done by adding the application
using `createGdcAppWithOwnStore`function. If the app is not using its store, then the `createGdcApp` function can be used.

```typescript
import { createGdcAppWithOwnStore } from "@gff/core";
import { AppContext, AppStore, id } from "@/features/projectsCenter/appApi";
import { ProjectsCenter } from "@/features/projectsCenter/ProjectsCenter";

export default createGdcAppWithOwnStore({
  App: ProjectsCenter,
  id: id,
  name: "Projects Center",
  version: "v1.0.0",
  requiredEntityTypes: [],
  store: AppStore,
  context: AppContext,
});

export const ProjectsCenterAppId: string = id;
```
The above code registers the application with the GDC Portal. The `createGdcAppWithOwnStore` function takes a number of arguments:

* `App`: React.ComponentType - the application component
* `id`: string - the id of the application
* `name`: string - the name of the application
* `version`: string - the version of the application
* `requiredEntityTypes`: string[] - the required entity types for the application
* `store`: Store - the store for the application
* `context`: Context - the context for the application

The required entity types are the entity types that the application requires to function. For example, the Mutation Frequency application
requires the `ssms` entity type. While this value is not currently used, it will be used in the future to determine if the application
can be used.

The other registration needed for your app is in
[packages/portal-proto/src/features/user-flow/workflow/registeredApps.tsx](https://github.com/NCI-GDC/gdc-frontend-framework/blob/f9ab9710450172978f5f588558cbdaa2d2301418/packages/portal-proto/src/features/user-flow/workflow/registeredApps.tsx)
This file contains an array of registered applications. For example the entry for the Project Center is:

```tsx
import ProjectsIcon from "public/user-flow/icons/crowd-of-users.svg";

...
{
  name: "Projects",
          icon: (
          <ProjectsIcon
                  width={64}
                  height={64}
                  viewBox="0 -20 128 128"
                  role="img"
                  aria-label="Projects icon"
          />
          ),
  tags: [],
          hasDemo: false,
          id: "Projects",
          countsField: "repositoryCaseCount",
          description:"View the Projects available within the GDC and select them for further exploration and analysis.",
},
...
```

The above code registers the Project Center application with the GDC Portal. The members of the object are:

*`name` is the name of the application
*`icon` is the icon as an SVH file, it size and position can be adjusted using the `width`, `height`, and `viewBox` properties
*`tags` are the tags for the application used for searching (which is not currently active)
*`hasDemo` is a boolean indicating if the application has a demo, if so the demo button will be shown
*`id` is the id of the application and needs to match the id of the application registered in the `createGdcAppWithOwnStore` function
*`countsField` is the field to use for the counts in the application, this is used to determine if the application can be used
*`description` is the description of the application
*`noDataTooltip` is the tooltip to show if the application has no data

When the app is registered, it will be available in the GDC Portal. The application can be accessed by clicking on the app card.
The visual elements of the card are:

![application_card.png](images%2Fapplication_card.png)

*Application Card and it's elements*

## Source code layout

While you are free to structure your application code as you with, the following is a recommended layout for your application's source code:

![source code layout](./images/app_source_code_layout_fig.png)

*Application source code layout*


# Appendix

## Using selectors and hooks

Although a complete guide to react hooks and selectors is out of the scope of this document, we will provide a brief overview
of how to use them for application development. For more information on hooks and selectors please see the
[React Hooks](https://reactjs.org/docs/hooks-intro.html). As we are using Redux-toolkit, we will be using the calls described in the [Redux Toolkit](https://redux-toolkit.js.org/tutorials/typescript) documentation.

### Selectors

Selectors are used to access the state of the GDC Portal's main redux store. Using selectors is the preferred method for
accessing the state of the GDC Portal. Selectors are functions that take the state as an argument and return a value.

```typescript
import {useCoreSelector,  selectCurrentCohort } from '@gff/core';

const currentCohort = useSelector(selectCurrentCohort);

```

The selector will return the current value of the item in the store. Consult the GDC V2 API documentation for a complete
list of selectors.

### Hooks

Fetching data from the GDC API is done via hooks. Hooks are functions that take arguments and return a value. The value
returned is typically a promise that resolves to the data requested. The GDC Portal provides a number of hooks for
fetching data from the GDC API. These hooks are located in the `@gff/core` package.

```typescript
import { useGeneSymbol } from '@gff/core';

const { data: geneSymbolDict, isSuccess } = useGeneSymbol(
        field === "genes.gene_id" ? facetValues.map((x) => x.toString()) : [],
);
```

GDC Portal hooks are design to work simlary to the RTL Query hooks. The hooks take arguments and return a object.
The object contains the data and the status of the query. The status of the query is stored in the `isSuccess` variable.
The data returned from the query is stored in the `data` variable. The object returned from a GDC hook is of the form:

```typescript
{
  data: any;
  isSuccess: boolean;
  isLoading: boolean;
  isError: boolean;
  error: Error;
}
```

where `data` is the data returned from the query, `isSuccess` is a boolean indicating if the query was successful, `isLoading`
is a boolean indicating if the query is currently loading, `isError` is a boolean indicating if the query resulted in an error,
and `error` is the error returned from the query.

## Querying the GDC API Directly

There may be cases where you need to query the GDC API directly. The GDC Portal provides a number of functions for querying
the GDC API. These functions are located in the `@gff/core` package. The functions are:
* `fetchGdcProjects` - fetches project data
* `fetchGdcAnnotations` - fetches annotation data
* `fetchGdcSsms` - fetches ssms data
* `fetchGdcCases` - fetches cases data
* `fetchGdcFiles` - fetches files data

which are wrappers around `fetchGdcEntities` function. The `fetchGdcEntities` function takes a number of arguments:

```typescript
export interface GdcApiRequest {
  readonly filters?: GqlOperation;
  readonly case_filters?: GqlOperation;
  readonly fields?: ReadonlyArray<string>;
  readonly expand?: ReadonlyArray<string>;
  readonly format?: "JSON" | "TSV" | "XML";
  readonly size?: number;
  readonly from?: number;
  readonly sortBy?: ReadonlyArray<SortBy>;
  readonly facets?: ReadonlyArray<string>;
}
```

There is also support for the GraphQL API. The `fetchGdcGraphQL` function takes two arguments:

```typescript
export const graphqlAPI = async <T>(
  query: string,
  variables: Record<string, unknown>,
): Promise<GraphQLApiResponse<T>> =>
```

where `query` is the GraphQL query and `variables` are the variables for the query.
