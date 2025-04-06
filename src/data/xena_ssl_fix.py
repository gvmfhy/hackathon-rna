#!/usr/bin/env python3
"""
SSL Certificate Fix for xenaPython

This module monkey-patches urllib.request to disable SSL verification
when used by xenaPython.

Usage:
    import xena_ssl_fix
    import xenaPython as xena
    # Now you can use xenaPython without SSL errors

Author: Hackathon-Bio Team
"""

import ssl
import urllib.request
import warnings

# Store the original urlopen function
original_urlopen = urllib.request.urlopen

def patch_urlopen():
    """
    Patch urllib.request.urlopen to disable SSL certificate verification.
    This is a workaround for SSL certificate verification errors.
    """
    def patched_urlopen(url, *args, **kwargs):
        # Create a context that doesn't verify certificates
        context = ssl._create_unverified_context()
        
        # Add the SSL context to the kwargs
        if 'context' not in kwargs:
            kwargs['context'] = context
        
        return original_urlopen(url, *args, **kwargs)
    
    # Apply the patch
    urllib.request.urlopen = patched_urlopen
    
    # Show warning
    warnings.warn(
        "SSL certificate verification has been disabled for xenaPython calls. "
        "This is a security risk and should be used only for development/testing.",
        UserWarning
    )

# Apply the patch immediately
patch_urlopen() 