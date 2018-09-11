# Data Security

To protect the privacy of research participants and support data integrity, the GDC requires user authorization and authentication for:

 - downloading controlled-access data
 - submitting data to the GDC

To perform these functions, GDC users must first obtain appropriate authorization via dbGaP and then authenticate via eRA Commons. The GDC sets user permissions at the project level according to dbGaP authorizations.

*See [Data Access Processes and Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools) to learn more about the difference between open-access and controlled-access data.*

## Authorization via dbGaP

Instructions for obtaining authorization via dbGaP are provided in [Obtaining Access to Controlled Data](https://gdc.cancer.gov/access-data/obtaining-access-controlled-data) and [Obtaining Access to Submit Data](https://gdc.cancer.gov/submit-data/obtaining-access-submit-data).

## Authentication via eRA Commons

The following authentication methods are supported by the GDC:

|GDC Tool|Authentication Method|
|----|----|
| GDC Data Portal | Log in using eRA Commons account|
| GDC Data Submission Portal | Log in using eRA Commons account |
| GDC Data Transfer Tool | Authentication Token |
| GDC API | Authentication Token |


### Authentication Tokens

The GDC Data Transfer Tool and the GDC API use tokens for authentication. GDC authentication tokens are alphanumeric strings of characters like this one:

	ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO


#### Obtaining A Token

Users can obtain authentication tokens from the [GDC Data Portal](https://portal.gdc.cancer.gov) and the [GDC Data Submission Portal](https://portal.gdc.cancer.gov/submission). See the [GDC Data Portal User's Guide](../../Data_Portal/Users_Guide/Authentication.md#gdc-authentication-tokens) and the [GDC Data Submission Portal User's Guide](../../Data_Submission_Portal/Users_Guide/Data_Submission_Process.md#authentication) for instructions.

#### Token Expiration

Tokens are valid for 30 days from the time of issue. Any request to the GDC API that uses an expired token will result in an error.

Tokens can be replaced at any time by downloading a new token, which will be valid for another 30 days.

## Checking User Permissions

Users can view the permissions granted to them by the GDC system as follows:

0. Log into the [GDC Data Portal](https://portal.gdc.cancer.gov) or the [GDC Data Submission Portal](https://portal.gdc.cancer.gov/submission) using your eRA Commons account.
0. Open the URL `https://portal.gdc.cancer.gov/auth/user` to see a JSON object that describes user permissions.
