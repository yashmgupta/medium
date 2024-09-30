# app.py

import streamlit as st
from data import research_summaries

# Initialize session state variables
if 'current_index' not in st.session_state:
    st.session_state.current_index = 0

if 'liked' not in st.session_state:
    st.session_state.liked = False

if 'query' not in st.session_state:
    st.session_state.query = ''

if 'filtered_summaries' not in st.session_state:
    st.session_state.filtered_summaries = research_summaries

# Search bar
st.title("📚 Research Reels")
query_input = st.text_input("Search topics...", value=st.session_state.query)

# Update the search query
if query_input != st.session_state.query:
    st.session_state.query = query_input
    # Filter the summaries based on the query
    st.session_state.filtered_summaries = [
        summary for summary in research_summaries
        if any(query_input.lower() in tag.lower() for tag in summary['tags'])
        or query_input.lower() in summary['title'].lower()
        or query_input.lower() in summary['summary'].lower()
    ]
    st.session_state.current_index = 0  # Reset index when search changes

# If no summaries found
if not st.session_state.filtered_summaries:
    st.write(f"No results found for '{st.session_state.query}'")
else:
    # Get the current summary
    current_summary = st.session_state.filtered_summaries[st.session_state.current_index]

    # Display the summary
    st.header(current_summary['title'])
    st.write(current_summary['summary'])
    st.write(f"**Tags:** {', '.join(current_summary['tags'])}")

    # Action buttons
    col1, col2, col3 = st.columns(3)

    with col1:
        if st.button('❤️ Like', key=f'like_{current_summary["id"]}'):
            st.session_state.liked = not st.session_state.liked

        if st.session_state.liked:
            st.write('You liked this!')

    with col2:
        comment = st.text_input("Add a comment...", key=f'comment_{current_summary["id"]}')
        if st.button("Submit Comment", key=f'submit_comment_{current_summary["id"]}'):
            st.write(f"Comment submitted: {comment}")

    with col3:
        if st.button('🔗 Share', key=f'share_{current_summary["id"]}'):
            st.write("Share functionality is not implemented yet.")

    # Navigation buttons
    nav1, nav2 = st.columns([1, 1])

    with nav1:
        if st.button('⬅️ Previous', key='prev'):
            if st.session_state.current_index > 0:
                st.session_state.current_index -= 1
                st.session_state.liked = False  # Reset liked status

    with nav2:
        if st.button('Next ➡️', key='next'):
            if st.session_state.current_index < len(st.session_state.filtered_summaries) - 1:
                st.session_state.current_index += 1
                st.session_state.liked = False  # Reset liked status

    # Pagination Indicator
    st.write(f"Summary {st.session_state.current_index + 1} of {len(st.session_state.filtered_summaries)}")
