# RAS Integration

## Coming Soon: Genomic Data Commons Integration with RAS

In the coming months, the GDC will be integrated with the NIH Researcher Auth Service (RAS). RAS provides a consistent NIH authentication and authorization framework for access to NIH open and controlled data. For controlled-access data, users will need an identity that meets the required Identity Assurance Level (IAL) before access can be granted ([Login.gov](https://login.gov/) or [ID.me](https://id.me/)). See below for information on how GDC users should prepare for these changes:

??? note "Preparing for RAS Integration"

    ### Preparing for RAS Integration

    Before the GDC moves to the new RAS login, we recommend preparing the identity you will use for controlled-access data. Completing this work beforehand will reduce the chance of an access interruption. GDC users should start at the Clinical and Translational Data Commons (CTDC) as they are currently integrated with RAS.  

    1. Go to the [CTDC Login Page](https://clinical.datacommons.cancer.gov/#/user/login) and choose the "Login with RAS" button. You'll be redirected to RAS and the appropriate identity provider.
    1. Use the Login.gov or ID.me links for GDC controlled-access data. These are the supported IAL2 identity providers.
    1. Complete the IAL2 identity proofing if your Login.gov or ID.me identity has not been verified at that level.
    1. Complete multi-factor authentication when prompted.
    1. For GDC controlled-access data, access permissions are managed through the Database of Genotypes and Phenotypes (dbGaP) and associated with your eRA Commons account. If RAS prompts you to link your eRA Commons identity, complete the identity-linking process.
    1. Review the RAS information-sharing consent screen. For non-NIH identities, you must grant the required consent for RAS to release identity and authorization information to the requesting application.
    1. Return to the CTDC application and confirm that the expected controlled-data authorization is available.

??? note "IAL2: what it means and why GDC users should prepare now"

    ### IAL2: what it means and why GDC users should prepare now

    IAL2 is an identity-proofing level defined in federal digital identity guidance. In practical terms, it means the identity used for controlled-data access has been verified against real-world evidence rather than relying only on self-asserted account information. RAS can detect whether the identity presented to a controlled-access repository meets the configured assurance requirement and can direct the user to complete additional proofing when needed.

    An eRA Commons username and password by itself is an IAL1 identity in RAS and is not an approved direct login method for a CADR requiring IAL2. Where needed, users can link an eRA identity to an IAL2 identity in RAS so that RAS can associate the user with relevant eRA/dbGaP information.

??? note "Choosing the IAL2 Identity"

    ### Choosing the IAL2 Identity

    | User Status | Recommended Starting Point |
    |---|---|
    | NIH Federal Staff| Login.gov or ID.me|
    | User with dbGaP permissions tied to eRA Commons| Use an IAL2 Login.gov or ID.me identity, then link eRA in RAS when prompted|
    | Extramural or university researcher | Login.gov or ID.me|
    | International researcher | Use Login.gov or ID.me, based on the identity-proofing options available for your circumstances|

??? note "Login.gov or ID.me Verification"

    ### Login.gov

    A user who needs to complete Login.gov identity verification should expect to provide identifying information and identity evidence. The current NIH RAS user guidance describes a process that can include a state-issued ID, a phone or computer camera, Social Security number, and a phone number or mailed verification route. Users who are already verified should use their existing account rather than creating a new one.

    ### ID.me

    A user who needs to complete ID.me identity verification should use an existing ID.me account if one already exists. ID.me permits only one verified account per individual. Verification may use a self-service path or, where available, a video-call path. Users should follow the RAS-to-ID.me flow so that the verified identity can be shared with NIH after consent.

    1. Use an existing ID.me account rather than creating a duplicate account.
    1. Set up more than one MFA method when practical to reduce the risk of losing access if one method is unavailable.
    1. If the ID.me flow displays HHS branding, that is expected for NIH use under the HHS portfolio.
    1. For ID.me account, identity-proofing, or document-upload problems, use [ID.me support](https://help.id.me/hc/en-us/articles/9198013708439-Verify-your-identity-with-IDMe-Self-Service). For access problems after successful verification, use the Cancer Research Data Commons (CRDC) partner support path. __Note:__ Once the GDC has fully transitioned to RAS, users should contact GDC support at: support@nci-gdc.datacommons.io

??? note "Linking an IAL2 identity to eRA Commons"

    Controlled-data authorizations for the GDC originate from dbGaP and are associated with an eRA Commons identity. When linking is required, RAS can link a Login.gov or ID.me identity to eRA so that authorization information associated with the linked identity can be represented in the RAS Passport. The user should perform this linking only through the RAS flow.

    1. Sign in to the CRDC partner system using your IAL2 identity. When RAS or the partner application indicates that eRA linking is required, continue to the Linked Identities page.
    1. Select the option to link eRA and sign in with the eRA Commons credentials associated with the relevant dbGaP authorization.
    1. Complete the process promptly. The RAS linking session is short and may require re-authentication if it expires.
    1. Confirm that the linked identity is shown, then return to the CRDC application.

### Frequently Asked Questions

??? note "What does IAL2 mean?"

    IAL2 is an identity-proofing level defined in federal digital identity guidance. In practical terms, it means the identity used for controlled-data access has been verified against real-world evidence rather than relying only on self-asserted account information. RAS can detect whether the identity presented to a controlled-access repository meets the configured assurance requirement and can direct the user to complete additional proofing when needed.

??? note "I have eRA Commons credentials, does that count as IAL2?"

    An eRA Commons username and password by itself is an IAL1 identity in RAS and is not an approved direct login method for a Controlled Access Data Repository (CADR) requiring IAL2. Where needed, users can link an eRA identity to an IAL2 identity in RAS so that RAS can associate the user with relevant eRA/dbGaP information.

??? note "Do I need to sign in through RAS to download open-access GDC data?"

    No, RAS authentication is not required to download open-access GDC data.

??? note "Basic Troubleshooting Steps"

    |Problem| Start Here|
    |---|---|
    |Cannot sign in to Login.gov or ID.me, duplicate account, or identity verification fails |Use the identity provider support channel. The GDC cannot help with problems from an external identity-provider account.|
    |Expected dbGaP authorization is missing |Confirm that the correct eRA identity is linked when linking is required and that the underlying dbGaP approval is active. Then contact partner support if the authorization is still absent.|
    |IAL2 verification succeeded, but the CRDC application still denies login| Contact the CRDC partner support team (support@nci-gdc.datacommons.io for the GDC) with the date/time, partner system, identity provider used, and any displayed error or transaction identifier. |
    |Open-access file cannot be retrieved |Report the file's UUID and error to support@nci-gdc.datacommons.io. Open-access retrieval should not depend on a dbGaP authorization Visa.|
