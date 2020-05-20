# GenomeFileUtil
This is the basic readme for this module. Include any usage or deployment instructions and links to other documentation here.

## Current Status

| Branch  | Build                                                              | Coverage                                                                         | LGTM Alerts                                                     |
| ------- | ------------------------------------------------------------------ | -------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| master  | [![Build Status](https://travis-ci.org/kbaseapps/GenomeFileUtil.svg?branch=master)](https://travis-ci.org/kbaseapps/GenomeFileUtil)  | [![Coverage Status](https://coveralls.io/repos/github/kbaseapps/GenomeFileUtil/badge.svg?branch=master)](https://coveralls.io/github/kbaseapps/GenomeFileUtil?branch=master)  | [![Total alerts](https://img.shields.io/lgtm/alerts/g/kbaseapps/GenomeFileUtil.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/kbaseapps/GenomeFileUtil/alerts/)  |


## Testing

The tests use data from KBase's CI environment. In `test_local/test.cfg`, set the following:

```
test_token={ci_dev_token}
kbase_endpoint=https://ci.kbase.us/services
auth_service_url=https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login
```
